"""
Fit Birch-Murnaghan equations of state to the bcc, fcc, and molecular
phases of solid chlorine from the FHI-aims EOS sweep.

Outputs:
    docs/eos_curves.png     E(V) per phase with BM fit lines
    docs/eos_data.csv       cleaned (V, E) data
    docs/eos_fits.csv       fitted parameters (V0, B0, B0') per phase

"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


RE_ENERGY = re.compile(
    r"Total energy of the DFT / Hartree-Fock s\.c\.f\. calculation\s*:\s*"
    r"([-\d.E+]+)\s*eV"
)
RE_VOLUME = re.compile(r"Unit cell volume\s*:\s*([\d.E+]+)\s*Angstrom\^3")
RE_NATOMS = re.compile(r"Number of atoms\s*:\s*(\d+)")


def parse_aims_out(path):
    text = Path(path).read_text(errors="replace")
    e = RE_ENERGY.search(text)
    v = RE_VOLUME.search(text)
    n = RE_NATOMS.search(text)
    return {
        "energy_eV": float(e.group(1)) if e else None,
        "volume_A3": float(v.group(1)) if v else None,
        "n_atoms":   int(n.group(1)) if n else None,
        "converged": "Self-consistency cycle converged" in text,
    }


def build_dataframe(repo_root):
    repo_root = Path(repo_root)
    rows = []
    for phase in ("bcc", "fcc", "cl"):
        eos_dir = repo_root / phase / "eos_m06l"
        if not eos_dir.exists():
            continue
        for a_dir in sorted(eos_dir.iterdir()):
            aims_out = a_dir / "scf" / "aims.out"
            if not aims_out.exists():
                continue
            try:
                lattice = float(a_dir.name.split("_", 1)[1])
            except (IndexError, ValueError):
                continue
            d = parse_aims_out(aims_out)
            d["phase"] = phase
            d["lattice_A"] = lattice
            rows.append(d)
    return pd.DataFrame(rows).sort_values(["phase", "lattice_A"]).reset_index(drop=True)


def filter_inconsistent(df, bimodal_gap_eV_per_atom=1.5):

    df = df.copy().reset_index(drop=True)
    df["E_per_atom"] = df["energy_eV"] / df["n_atoms"]
    df["V_per_atom"] = df["volume_A3"] / df["n_atoms"]
    keep_mask = np.zeros(len(df), dtype=bool)

    for phase, group in df.groupby("phase"):
        idx = group.index.to_numpy()
        E = group["E_per_atom"].to_numpy()
        if len(E) < 3:
            keep_mask[idx] = True
            continue
        # Sort by E_per_atom
        order = np.argsort(E)
        sorted_E = E[order]
        gaps = np.diff(sorted_E)
        max_gap_pos = int(np.argmax(gaps))
        max_gap = float(gaps[max_gap_pos])

        if max_gap > bimodal_gap_eV_per_atom:
            # keep only the lower  cluster
            cutoff = sorted_E[max_gap_pos] + max_gap / 2.0
            keep_mask[idx] = E <= cutoff
        else:
            keep_mask[idx] = True

    kept = df[keep_mask].copy()
    dropped = df[~keep_mask].copy()
    return kept, dropped

# Birch-Murnaghan

EV_PER_A3_TO_GPA = 160.21766208


def birch_murnaghan_E(V, E0, V0, B0, B0p):
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + (9.0 * V0 * B0 / 16.0) * (
        (eta - 1.0) ** 3 * B0p + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )


def fit_bm(V, E):
    V = np.asarray(V, float)
    E = np.asarray(E, float)
    a, b, c = np.polyfit(V, E, 2)
    V0_guess = -b / (2 * a)
    p0 = (c - b * b / (4 * a), V0_guess, max(2 * a * V0_guess, 1e-3), 4.0)
    popt, _ = curve_fit(birch_murnaghan_E, V, E, p0=p0, maxfev=20_000)
    rmse = float(np.sqrt(np.mean((E - birch_murnaghan_E(V, *popt)) ** 2)))
    return (*popt, rmse)


# Plot
PHASE_STYLE = {
    "cl":  {"label": r"Cl$_2$", "color": "#ffa600", "marker": "o"},
    "fcc": {"label": "FCC",     "color": "#00d36b", "marker": "s"},
    "bcc": {"label": "BCC",     "color": "#e63374", "marker": "^"},
}


def plot_eos(df_kept, fits, out_path):
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    for phase in ("cl", "fcc", "bcc"):
        if phase not in fits:
            continue
        g = df_kept[df_kept["phase"] == phase].sort_values("V_per_atom")
        s = PHASE_STYLE[phase]
        # Data points
        ax.plot(g["V_per_atom"], g["E_per_atom"],
                marker=s["marker"], linestyle="none", markersize=7,
                color=s["color"], label=s["label"])
        # BM fit curve (per-atom)
        E0, V0_cell, B0, B0p, _ = fits[phase]
        n_atoms = int(g["n_atoms"].iloc[0])
        V_grid_cell = np.linspace(g["volume_A3"].min(), g["volume_A3"].max(), 300)
        E_grid_cell = birch_murnaghan_E(V_grid_cell, E0, V0_cell, B0, B0p)
        ax.plot(V_grid_cell / n_atoms, E_grid_cell / n_atoms,
                color=s["color"], linewidth=2, alpha=0.7)
    ax.set_xlabel(r"Volume per atom [Å$^3$]", fontsize=14)
    ax.set_ylabel("Energy per atom [eV]", fontsize=14)
    ax.set_title("Birch-Murnaghan equation of state, solid Cl", fontsize=14)
    ax.legend(fontsize=12, frameon=True)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)


def main():
    repo_root = Path(__file__).resolve().parent
    docs = repo_root / "docs"
    docs.mkdir(exist_ok=True)
 
    df = build_dataframe(repo_root)
    if df.empty:
        raise SystemExit("No aims.out files found. Run from the repo root.")
 
    df = df.dropna(subset=["energy_eV", "volume_A3", "n_atoms"])
    df = df[df["converged"]].copy()
    df["n_atoms"] = df["n_atoms"].astype(int)
 
    kept, _ = filter_inconsistent(df)
 
    fits = {}
    fit_rows = []
    for phase, group in kept.groupby("phase"):
        if len(group) < 4:
            continue
        E0, V0, B0, B0p, rmse = fit_bm(group["volume_A3"], group["energy_eV"])
        n_atoms = int(group["n_atoms"].iloc[0])
        fits[phase] = (E0, V0, B0, B0p, rmse)
        fit_rows.append({
            "phase": phase,
            "n_atoms_per_cell": n_atoms,
            "V0_per_atom_A3": V0 / n_atoms,
            "B0_GPa": B0 * EV_PER_A3_TO_GPA,
            "B0_prime": B0p,
            "rmse_meV": rmse * 1000,
        })
 
    kept.to_csv(docs / "eos_data.csv", index=False)
    pd.DataFrame(fit_rows).to_csv(docs / "eos_fits.csv", index=False)
    plot_eos(kept, fits, docs / "eos_curves.png")
 
 
if __name__ == "__main__":
    main()