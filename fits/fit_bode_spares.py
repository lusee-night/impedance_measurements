#!/usr/bin/env python3
"""
Re-analysis of the Bode 500 sweeps of the four SPARE preamps (fmpre0/2/5/7),
BODE500_BNL/data/<preamp>_impedance_bode*.csv, 10 sweeps each.

The fitting formula is taken as given from Toki's presentation
(BODE500_BNL/2026_06_01_LuSEENight_JFETAnalysis.pdf):

    Re Z = 1/(w^2 Rp C^2) + Rs + a w^n
    Im Z = w L - 1/(w C)

fitted over 1-100 MHz to the mean of the 10 sweeps, weighted by the standard
error of the mean.  Three fits are done per preamp:

  * "Im-only":  C, L               from Im Z
  * "Re-only":  Rs, n, a, Rp       from Re Z  (C fixed from the Im-only fit)
  * "joint":    all six parameters from Re Z and Im Z together

and the results are compared with the table in the presentation.  Two
weightings are run: "SE" (1/SE^2 of the sweep mean) and "unweighted"; the
latter reproduces the presentation's C, L (and their errors) closely, so that
is presumably what was used there.  Rs, n, a and Rp are strongly degenerate
(Rs + a w^n and 1/(w^2 Rp C^2) trade off against each other), so compare
R_L(100 MHz) = a (2 pi 1e8)^n rather than a and n separately.

Run:  python3 fits/fit_bode_spares.py
Writes fits/bode_spares_refit.png.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BODE_DIR = os.path.join(ROOT, "BODE500_BNL", "data")
PREAMPS = ["fmpre0", "fmpre2", "fmpre5", "fmpre7"]
FMIN, FMAX = 1e6, 100e6

# Table "Fitted Parameters - All 4 Preamps" from the presentation (values, errors).
PRESENTATION = {
    "fmpre0": {"C_pf": (33.79, 0.03), "L_nh": (35.0, 1.8), "Rs_ohm": (0.21, 0.01), "n": (1.39, 0.00),
               "a": (5.2e-12, None), "Rp_kohm": (1101, 18)},
    "fmpre2": {"C_pf": (35.32, 0.04), "L_nh": (36.4, 2.0), "Rs_ohm": (2.16, 0.00), "n": (2.26, 0.01),
               "a": (6.9e-20, None), "Rp_kohm": (331, 2)},
    "fmpre5": {"C_pf": (34.89, 0.04), "L_nh": (34.7, 1.9), "Rs_ohm": (3.71, 0.00), "n": (4.35, 0.01),
               "a": (2.2e-38, None), "Rp_kohm": (288, 1)},
    "fmpre7": {"C_pf": (36.85, 0.05), "L_nh": (38.5, 2.5), "Rs_ohm": (1.65, 0.01), "n": (1.88, 0.01),
               "a": (1.9e-16, None), "Rp_kohm": (375, 2)},
}
COL = {"fmpre0": "#e87ba4", "fmpre2": "#008300", "fmpre5": "#6250d6", "fmpre7": "#e34948"}


def load_bode(preamp):
    files = sorted(glob.glob(os.path.join(BODE_DIR, f"{preamp}_impedance_bode*.csv")))
    A = np.array([pd.read_csv(f, sep="\t", encoding="utf-8-sig").iloc[:, :3].values for f in files])
    k = len(files)
    return dict(f=A[0, :, 0], re=A[:, :, 1].mean(0), re_se=A[:, :, 1].std(0, ddof=1) / np.sqrt(k),
                im=A[:, :, 2].mean(0), im_se=A[:, :, 2].std(0, ddof=1) / np.sqrt(k), nsweeps=k)


# parameter vector p = [C_pf, L_nh, Rs_ohm, n, log10(a), Rp_kohm]
def re_model(f, C_pf, Rs_ohm, n, log10a, Rp_kohm):
    w = 2 * np.pi * f
    C = C_pf * 1e-12
    return 1.0 / (w**2 * Rp_kohm * 1e3 * C**2) + Rs_ohm + 10.0**log10a * w**n


def im_model(f, C_pf, L_nh):
    w = 2 * np.pi * f
    return w * L_nh * 1e-9 - 1.0 / (w * C_pf * 1e-12)


def _fit(resid_fn, p0, bounds):
    r = least_squares(resid_fn, p0, bounds=bounds, x_scale="jac", max_nfev=20000)
    chi2 = float(np.sum(r.fun**2))
    dof = len(r.fun) - len(p0)
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J)
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(len(p0), np.nan)
    return r.x, err, chi2, dof


def fit_preamp(b, weighting="SE"):
    sel = (b["f"] >= FMIN) & (b["f"] <= FMAX)
    f, re, re_se, im, im_se = (b[k][sel] for k in ("f", "re", "re_se", "im", "im_se"))
    if weighting == "unweighted":
        re_se, im_se = np.ones_like(re), np.ones_like(im)

    # Im-only: C, L
    pI, eI, chi2I, dofI = _fit(lambda p: (im_model(f, *p) - im) / im_se, [35.0, 35.0],
                               ([1, -1e3], [100, 1e3]))
    # Re-only: Rs, n, log10a, Rp_kohm with C fixed
    C_fix = pI[0]
    pR, eR, chi2R, dofR = _fit(lambda p: (re_model(f, C_fix, *p) - re) / re_se, [1.0, 2.0, -16.0, 400.0],
                               ([0, 0.1, -80, 1], [100, 10, 0, 1e5]))
    # joint
    def rj(p):
        C, L, Rs, n, la, Rp = p
        return np.concatenate([(re_model(f, C, Rs, n, la, Rp) - re) / re_se,
                               (im_model(f, C, L) - im) / im_se])
    p0 = [pI[0], pI[1], pR[0], pR[1], pR[2], pR[3]]
    pJ, eJ, chi2J, dofJ = _fit(rj, p0, ([1, -1e3, 0, 0.1, -80, 1], [100, 1e3, 100, 10, 0, 1e5]))
    if weighting == "unweighted":  # scale errors by the residual variance
        for (pp, ee, c2, dd) in ((pI, eI, chi2I, dofI), (pR, eR, chi2R, dofR), (pJ, eJ, chi2J, dofJ)):
            ee *= np.sqrt(c2 / dd)
    return dict(Im=(pI, eI, chi2I, dofI), Re=(pR, eR, chi2R, dofR), joint=(pJ, eJ, chi2J, dofJ), f=f)


def fmt(v, e, unit=""):
    return f"{v:.2f} ± {e:.2f}{unit}" if e is not None and np.isfinite(e) else f"{v:.2f}{unit}"


def main():
    data = {p: load_bode(p) for p in PREAMPS}
    allfits = {w: {p: fit_preamp(data[p], w) for p in PREAMPS} for w in ("SE", "unweighted")}
    for weighting, fits in allfits.items():
        print_table(fits, weighting)
    plot(data, allfits["unweighted"])
    return allfits


def print_table(fits, weighting):
    print(f"\n{'='*80}\nFit range {FMIN/1e6:.0f}-{FMAX/1e6:.0f} MHz, mean of 10 sweeps, weighting: {weighting}\n")
    hdr = f"{'preamp':7s} {'quantity':9s} {'presentation':>18s} {'Im-only / Re-only':>20s} {'joint':>20s}"
    for p in PREAMPS:
        F = fits[p]
        pI, eI, chi2I, dofI = F["Im"]
        pR, eR, chi2R, dofR = F["Re"]
        pJ, eJ, chi2J, dofJ = F["joint"]
        P = PRESENTATION[p]
        print(hdr)
        rows = [
            ("C (pF)", P["C_pf"], (pI[0], eI[0]), (pJ[0], eJ[0])),
            ("L (nH)", P["L_nh"], (pI[1], eI[1]), (pJ[1], eJ[1])),
            ("Rs (Ohm)", P["Rs_ohm"], (pR[0], eR[0]), (pJ[2], eJ[2])),
            ("n", P["n"], (pR[1], eR[1]), (pJ[3], eJ[3])),
            ("Rp (kOhm)", P["Rp_kohm"], (pR[3], eR[3]), (pJ[5], eJ[5])),
        ]
        for name, pres, sep, joint in rows:
            pv = f"{pres[0]:.2f} ± {pres[1]:.2f}" if pres[1] is not None else f"{pres[0]:.3g}"
            print(f"{p:7s} {name:9s} {pv:>18s} {fmt(*sep):>20s} {fmt(*joint):>20s}")
        a_pres = P["a"][0]
        print(f"{p:7s} {'a (Ohm s^n)':9s} {a_pres:>18.1e} {10**pR[2]:>20.1e} {10**pJ[4]:>20.1e}")
        # the (a, n) pair is degenerate; compare the resulting R_L at 100 MHz instead
        w100 = 2 * np.pi * 100e6
        print(f"{p:7s} {'R_L@100MHz':9s} {a_pres * w100**P['n'][0]:>17.2f}Ω {10**pR[2] * w100**pR[1]:>19.2f}Ω "
              f"{10**pJ[4] * w100**pJ[3]:>19.2f}Ω")
        print(f"{p:7s} chi2/dof: Im-only {chi2I:.0f}/{dofI}, Re-only {chi2R:.0f}/{dofR}, joint {chi2J:.0f}/{dofJ}\n")


def plot(data, fits):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(17, 7), sharex=True)
    for j, p in enumerate(PREAMPS):
        b, F = data[p], fits[p]
        f = F["f"]
        sel = (b["f"] >= FMIN) & (b["f"] <= FMAX)
        pJ = F["joint"][0]
        ff = np.logspace(np.log10(FMIN), np.log10(FMAX), 400)
        P = PRESENTATION[p]
        ax = axes[0, j]
        ax.errorbar(f, b["re"][sel], yerr=b["re_se"][sel], fmt="o", ms=3, color="#52514e", lw=0.8, label="data (mean ± SE)")
        ax.plot(ff, re_model(ff, pJ[0], pJ[2], pJ[3], pJ[4], pJ[5]), color=COL[p], lw=2, label="joint refit")
        ax.plot(ff, re_model(ff, P["C_pf"][0], P["Rs_ohm"][0], P["n"][0], np.log10(P["a"][0]), P["Rp_kohm"][0]),
                color="#0b0b0b", lw=1.2, ls="--", label="presentation params")
        ax.set_title(p)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(FMIN, FMAX)
        ax.set_ylim(0.3, 30)
        ax.grid(True, which="both", color="#dddddd", lw=0.6)
        if j == 0:
            ax.set_ylabel("Re Z [Ω]")
            ax.legend(fontsize=8)
        ax = axes[1, j]
        ax.errorbar(f, -b["im"][sel], yerr=b["im_se"][sel], fmt="o", ms=3, color="#52514e", lw=0.8)
        ax.plot(ff, -im_model(ff, pJ[0], pJ[1]), color=COL[p], lw=2)
        ax.plot(ff, -im_model(ff, P["C_pf"][0], P["L_nh"][0]), color="#0b0b0b", lw=1.2, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("frequency [Hz]")
        ax.grid(True, which="both", color="#dddddd", lw=0.6)
        if j == 0:
            ax.set_ylabel("−Im Z [Ω]")
    fig.suptitle("Bode 500 spare preamps, 1–100 MHz: data, joint unweighted refit (solid), and presentation parameters (dashed)")
    fig.tight_layout()
    out = os.path.join(HERE, "bode_spares_refit.png")
    fig.savefig(out, dpi=120)
    print(f"Figure written to {out}")


if __name__ == "__main__":
    main()
