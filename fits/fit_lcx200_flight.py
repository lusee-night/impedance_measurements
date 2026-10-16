#!/usr/bin/env python3
"""
Fit the input capacitance of the four FLIGHT preamps from the LCX200 LCR-meter
data taken at SDL (LCX200_SSL_SDL/data/ZJFET_250930.csv), and compare with the
Bode 500 sweeps of the four SPARE preamps taken at BNL (BODE500_BNL/data/).

Model (same as in preamp_fits.py, from Toki's analysis):
    Re Z = 1/(w^2 Rp C^2) + Rs + a w^n
    Im Z = w L - 1/(w C)
For the flight units only C is fitted; L, Rs, a, n, Rp are fixed to the
spare-preamp averages (see spare_average()).

Run:  python3 fits/fit_lcx200_flight.py
Writes fits/lcx200_flight_capacitance.png and prints the FLIGHT_PREAMP_FITS block.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from preamp_fits import SPARE_PREAMP_FITS  # noqa: E402

# Flight model preamp id -> ADC channel (README.md, "Flight model connections").
# The LCX200 CSV "Preamp" column (0..3) is taken to be the ADC channel.
CHANNEL_TO_FLIGHT_PREAMP = {0: "fmpre6", 1: "fmpre3", 2: "fmpre4", 3: "fmpre1"}

LCX_CSV = os.path.join(ROOT, "LCX200_SSL_SDL", "data", "ZJFET_250930.csv")
BODE_DIR = os.path.join(ROOT, "BODE500_BNL", "data")

# Fixed categorical colours, one per entity (never re-used).
COL = {"fmpre6": "#2a78d6", "fmpre3": "#eb6834", "fmpre4": "#1baf7a", "fmpre1": "#eda100",
       "fmpre0": "#e87ba4", "fmpre2": "#008300", "fmpre5": "#6250d6", "fmpre7": "#e34948"}


# ----------------------------------------------------------------------------
def model_Z(f, C_pf, L_nh, Rs_ohm, a, n, Rp_ohm):
    w = 2 * np.pi * f
    C = C_pf * 1e-12
    re = 1.0 / (w**2 * Rp_ohm * C**2) + Rs_ohm + a * w**n
    im = w * L_nh * 1e-9 - 1.0 / (w * C)
    return re + 1j * im


def spare_average(fits=SPARE_PREAMP_FITS, fgrid=np.logspace(6, 8, 200)):
    """Average the non-C parameters over the spare preamps.

    L, Rs: arithmetic mean.  Rp: harmonic mean (its Re Z contribution is
    proportional to 1/Rp, so this averages the actual contribution).
    (a, n): the individual power laws a*w^n have wildly different (a, n) pairs
    that cannot be averaged term by term; instead average the curves
    R_L(w) = a*w^n over 1-100 MHz and refit a single power law.
    """
    ps = list(fits.values())
    L = np.mean([p["L_nh"] for p in ps])
    Rs = np.mean([p["Rs_ohm"] for p in ps])
    Rp = 1.0 / np.mean([1.0 / p["Rp_ohm"] for p in ps])
    w = 2 * np.pi * fgrid
    RL = np.mean([p["a"] * w ** p["n"] for p in ps], axis=0)
    n, loga = np.polyfit(np.log(w), np.log(RL), 1)
    return {"L_nh": float(L), "Rs_ohm": float(Rs), "n": float(n),
            "a": float(np.exp(loga)), "Rp_ohm": float(Rp)}, RL, fgrid


# ----------------------------------------------------------------------------
def load_lcx200():
    d = pd.read_csv(LCX_CSV, encoding="utf-8-sig")
    d["Z_ohm"] = d["Z"] * 1e3          # kOhm -> Ohm
    d["f_hz"] = d["Frequency"] * 1e3   # kHz -> Hz
    d["Re"] = d["Z_ohm"] * np.cos(np.radians(d["Phase"]))
    d["Im"] = d["Z_ohm"] * np.sin(np.radians(d["Phase"]))
    g = d.groupby(["Preamp", "f_hz"])
    s = g.agg(n=("Z_ohm", "size"),
              Zmag=("Z_ohm", "mean"), Zmag_se=("Z_ohm", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
              phase=("Phase", "mean"), phase_se=("Phase", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
              Im=("Im", "mean"), Im_se=("Im", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
              Re=("Re", "mean"), Re_se=("Re", lambda x: x.std(ddof=1) / np.sqrt(len(x)))).reset_index()
    s["w"] = 2 * np.pi * s["f_hz"]
    s["C_app_pf"] = 1e12 / (s["w"] * s["Zmag"])
    s["C_app_se"] = s["C_app_pf"] * s["Zmag_se"] / s["Zmag"]
    return s


def load_bode(preamp):
    files = sorted(glob.glob(os.path.join(BODE_DIR, f"{preamp}_impedance_bode*.csv")))
    arrs = [pd.read_csv(f, sep="\t", encoding="utf-8-sig").iloc[:, :3].values for f in files]
    A = np.array(arrs)
    f = A[0, :, 0]
    re, im = A[:, :, 1], A[:, :, 2]
    zmag = np.hypot(re, im)
    k = len(files)
    return dict(f=f, re=re.mean(0), re_se=re.std(0, ddof=1) / np.sqrt(k),
                im=im.mean(0), im_se=im.std(0, ddof=1) / np.sqrt(k),
                zmag=zmag.mean(0), zmag_se=zmag.std(0, ddof=1) / np.sqrt(k), nsweeps=k)


# ----------------------------------------------------------------------------
def fit_C(sub, fixed):
    """Weighted LSQ of |Z| with C the only free parameter.  Returns
    (C, stat err, err scaled by sqrt(chi2/dof), chi2, dof)."""
    f, y, sy = sub["f_hz"].values, sub["Zmag"].values, sub["Zmag_se"].values

    def fn(f, C_pf):
        return np.abs(model_Z(f, C_pf, **fixed))

    if len(f) == 1:  # one point: solve directly, propagate error
        popt, _ = curve_fit(fn, f, y, p0=[35.0])
        C = popt[0]
        dCdZ = -C / y[0]  # C ~ 1/(w|Z|)
        err = abs(dCdZ * sy[0])
        return C, err, err, 0.0, 0
    popt, pcov = curve_fit(fn, f, y, p0=[35.0], sigma=sy, absolute_sigma=True)
    chi2 = float(np.sum(((y - fn(f, *popt)) / sy) ** 2))
    dof = len(f) - 1
    err = float(np.sqrt(pcov[0, 0]))
    return float(popt[0]), err, err * np.sqrt(max(chi2 / dof, 1.0)), chi2, dof


def main():
    fixed, RL_avg, fgrid = spare_average()
    print("Spare-preamp averaged parameters (non-C):")
    for k, v in fixed.items():
        print(f"   {k:7s} = {v:.4g}")
    print("   (fmpre0/2/5/7 L_nh: ", [p["L_nh"] for p in SPARE_PREAMP_FITS.values()], ")")

    s = load_lcx200()
    print("\nLCX200 flight data, per channel and frequency (5 repeats each):")
    print(s[["Preamp", "f_hz", "n", "Zmag", "Zmag_se", "phase", "phase_se", "C_app_pf", "C_app_se"]]
          .to_string(index=False, float_format=lambda x: f"{x:10.3f}"))

    selections = {"all 5 freqs": lambda df: df,
                  ">=320 kHz": lambda df: df[df.f_hz >= 320e3],
                  "1 MHz only": lambda df: df[df.f_hz >= 1e6]}
    results = {}
    print("\nC fits (|Z|, full model with spare-averaged L,Rs,a,n,Rp fixed):")
    print(f"{'ch':>2s} {'preamp':>7s} | " + " | ".join(f"{k:^34s}" for k in selections))
    for ch in sorted(s.Preamp.unique()):
        pre = CHANNEL_TO_FLIGHT_PREAMP[ch]
        sub_all = s[s.Preamp == ch]
        row = []
        for key, selfn in selections.items():
            C, e, es, chi2, dof = fit_C(selfn(sub_all), fixed)
            results[(pre, key)] = (C, e, es, chi2, dof)
            row.append(f"C={C:6.2f} ±{e:4.2f} (±{es:4.2f} scaled) chi2/dof={chi2:7.1f}/{dof}")
        print(f"{ch:2d} {pre:>7s} | " + " | ".join(row))

    # ------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax_c, ax_ph, ax_re = axes
    for ax in axes:
        ax.set_xscale("log")
        ax.grid(True, which="both", color="#dddddd", lw=0.6)
        ax.set_xlabel("frequency [Hz]")
        ax.set_xlim(8e4, 1.3e6)
    # spares from Bode
    for pre in SPARE_PREAMP_FITS:
        b = load_bode(pre)
        sel = (b["f"] >= 9e4) & (b["f"] <= 1.1e6)
        f = b["f"][sel]
        capp = 1e12 / (2 * np.pi * f * b["zmag"][sel])
        ax_c.plot(f, capp, color=COL[pre], lw=1.5, ls="--", label=f"{pre} spare (Bode)")
        ax_ph.plot(f, np.degrees(np.arctan2(b["im"][sel], b["re"][sel])), color=COL[pre], lw=1.5, ls="--")
        ax_re.plot(f, b["re"][sel], color=COL[pre], lw=1.5, ls="--")
    # flight from LCX200
    for ch in sorted(s.Preamp.unique()):
        pre = CHANNEL_TO_FLIGHT_PREAMP[ch]
        sub = s[s.Preamp == ch]
        ax_c.errorbar(sub.f_hz, sub.C_app_pf, yerr=sub.C_app_se, color=COL[pre], marker="o", ms=6,
                      lw=1.5, capsize=3, label=f"{pre} flight CH{ch} (LCX200)")
        ax_ph.errorbar(sub.f_hz, sub.phase, yerr=sub.phase_se, color=COL[pre], marker="o", ms=6, lw=1.5, capsize=3)
        ax_re.errorbar(sub.f_hz, sub.Re, yerr=sub.Re_se, color=COL[pre], marker="o", ms=6, lw=1.5, capsize=3)
    ax_c.set_ylabel("apparent C = 1/(ω|Z|)  [pF]")
    ax_c.set_title("Apparent input capacitance")
    ax_ph.set_ylabel("phase of Z [deg]")
    ax_ph.set_title("Phase")
    ax_ph.axhline(-90, color="#888888", lw=0.8)
    ax_re.set_ylabel("Re Z [Ω]")
    ax_re.set_title("Real part")
    ax_re.axhline(0, color="#888888", lw=0.8)
    ax_c.legend(fontsize=8, ncol=1, loc="upper right")
    fig.suptitle("LuSEE-Night preamp input impedance: flight units (LCX200 @ SDL, points) vs spares (Bode 500 @ BNL, dashed)")
    fig.tight_layout()
    out = os.path.join(HERE, "lcx200_flight_capacitance.png")
    fig.savefig(out, dpi=130)
    print(f"\nFigure written to {out}")

    # ------------------------------------------------------------- dict block
    print("\nFLIGHT_PREAMP_FITS block (C from the 1 MHz point; see fits/README notes):")
    print("FLIGHT_PREAMP_FITS = {")
    for ch in sorted(s.Preamp.unique()):
        pre = CHANNEL_TO_FLIGHT_PREAMP[ch]
        C = results[(pre, "1 MHz only")][0]
        print(f'    "{pre}": {{"C_pf": {C:.2f}, "L_nh": {fixed["L_nh"]:.1f}, "Rs_ohm": {fixed["Rs_ohm"]:.2f}, '
              f'"n": {fixed["n"]:.2f},\n               "a": {fixed["a"]:.1e}, "Rp_ohm": {fixed["Rp_ohm"]/1e3:.0f}e3}},  # CH{ch}')
    print("}")
    return results, fixed


if __name__ == "__main__":
    main()
