# Preamplifier input impedance fits

There were two sets of measurements taken to determine the LuSEE-Night impedances.

First, some measurements have been taken at SDL using the actual flight hardware in final configuration, including the real flight hardware, actual cable length, etc. But these only go 1MHz. See [LCX200_SSL_SDL/](LCX200_SSL_SDL/).

Second, more measurements have been taken at BNL using BODE and going all the way to 50MHz, but on the flight spares. Toki then analyzed this data to find a fitting function that fits well to 50MHz. See [BODE500_BNL/](BODE500_BNL/).

Strategy: we will use C_pf measurements from LCX200_SSL_SDL for real flight hardware and the average values of other parameters (they vary very little anyway and are perturbative correction) as our impedance model.


# Preamp mapping

## Flight model connections

We list the flight model preamp id -> spectrometer connector -> ADC channel

fmpre6 -> J4 -> CH0
fmpre3 -> J5 -> CH1
fmpre4 -> J6 -> CH2
fmpre1 -> J7 -> CH3

## BNL flat sat preamps

These are the preamps that are connected to the BNL flat sat: fmpre0, fmpre2, fmpre5, fmpre7 and for which we BODE impendace measurement.




# Fits directory (`fits/`)

`fits/preamp_fits.py` holds the impedance model parameters used downstream:

    Re Z = 1/(w^2 Rp C^2) + Rs + a w^n,     Im Z = w L - 1/(w C)

with two dictionaries, `SPARE_PREAMP_FITS` (fmpre0/2/5/7, from the Bode data) and
`FLIGHT_PREAMP_FITS` (fmpre6/3/4/1, from the LCX200 data), plus
`FLIGHT_CHANNEL_TO_PREAMP` (CH0..CH3 -> fmpre id, from the table above).
The LCX200 CSV labels preamps 0..3; these are taken to be the ADC channels CH0..CH3.

## Spare preamps: check against Toki's presentation

`SPARE_PREAMP_FITS` was checked against the table "Fitted Parameters - All 4
Preamps" in `BODE500_BNL/2026_06_01_LuSEENight_JFETAnalysis.pdf`. One typo was
found and fixed: fmpre7 had L = 8.5 nH instead of 38.5 nH. Everything else matches.

`fits/fit_bode_spares.py` refits the Bode data (mean of the 10 sweeps per preamp,
1-100 MHz) with the formula above taken as given, and prints a comparison with the
presentation. Findings:

- An unweighted least-squares fit reproduces the presentation's C to 0.05 pF and L
  to about 1 nH, with the same 1-2 nH errors on L, so that is presumably what was
  used. A fit weighted by the standard error of the 10 sweeps moves C by up to
  0.2 pF and L by up to 2 nH (the model is not perfect at the level of the very
  small sweep-to-sweep scatter, chi^2/dof is in the thousands).
- The Re Z parameters (Rs, a, n, Rp) are strongly degenerate with each other and
  depend on the weighting; they should not be compared one by one. The physically
  relevant combination, the inductor loss R_L(100 MHz) = a (2 pi 1e8)^n, agrees
  with the presentation to about 1 Ohm. For fmpre0 the Re Z data below 5 MHz are
  very noisy and Rp is essentially unconstrained (the presentation's 1101 kOhm
  should not be taken literally). None of this matters in practice: Re Z is a few
  Ohm against a reactance of tens of Ohm to several kOhm.
- Figure: `fits/bode_spares_refit.png`.

## Flight preamps: analysis of the LCX200 data

`fits/fit_lcx200_flight.py` analyses `LCX200_SSL_SDL/data/ZJFET_250930.csv`
(5 repeat measurements per preamp at 100, 180, 320, 560 and 1000 kHz; |Z| in kOhm
and phase in degrees). For each preamp and frequency it takes the mean of the 5
repeats; the quoted errors are the standard error of that mean (scatter / sqrt(5)).
It then fits C with the full model, with L, Rs, a, n, Rp fixed to the spare
averages (see below), and writes `fits/lcx200_flight_capacitance.png`, which
overlays the flight data with the Bode data of the spares over the same band.

The flight data do not look like the spares over 100 kHz - 1 MHz:

- The spares (Bode, calibrated open/short/load) are a flat 34-37 pF capacitor
  with phase -90 +- 1 deg over the whole band.
- The LCX200 flight data have a phase of about -100 deg at every frequency,
  i.e. a negative real part of several kOhm, which is unphysical for a passive
  input and inconsistent with the -89.4 deg the Bode sees at 1 MHz on the spares.
- The apparent capacitance C = 1/(w |Z|) rises steadily from about 33 pF at
  1 MHz to 36 pF at 560 kHz, 38 pF at 320 kHz, 40-42 pF at 180 kHz and 44-47 pF
  at 100 kHz, the same for all four channels. This is the frequency-dependent
  capacitance Kaja saw. Fitting a pure capacitor (or the full model) to all five
  frequencies gives chi^2/dof > 400, which is Emi's "not statistically
  acceptable" result.
- The procedure (LST-SPT-0110) sets the LCX200 to its manual 3 kOhm range and
  says to keep pressing Trig until the RANGE warning goes away. Only the 1 MHz
  point (|Z| = 4.8 kOhm) is near that range; 100 kHz (36 kOhm) is 12x out of
  range, and the repeat scatter there is 5%, versus 0.1-0.3% at 1 MHz.

Conclusion: the frequency dependence and the phase offset are taken to be an
instrument effect of the out-of-range LCR measurement, not a property of the
flight preamps (same design as the spares, which show no such behaviour).
Therefore `C_pf` in `FLIGHT_PREAMP_FITS` is taken from the 1 MHz point only,
essentially C = 1/(w |Z|) at 1 MHz:

    CH0 fmpre6  33.12 +- 0.02 pF
    CH1 fmpre3  33.29 +- 0.11 pF
    CH2 fmpre4  33.44 +- 0.07 pF
    CH3 fmpre1  33.17 +- 0.06 pF

These sit at the low end of the spare range (33.8-36.9 pF). The statistical
errors (from the 5 repeats) are tiny; a realistic systematic error is a couple
of pF, e.g. the 560 kHz point is 2.7 pF higher. For reference, fitting all five
frequencies instead gives 33.4 / 37.7 / 35.6 / 35.3 pF for CH0..CH3, and Emi's
earlier recommendation for ch0 was 39 +- 2 pF. Switching to the all-frequency fit
is a one-line change in the script (the `"1 MHz only"` selection in the printed
block).

The other parameters of `FLIGHT_PREAMP_FITS` are averages over the four spares:
L and Rs arithmetic means, Rp the harmonic mean (which averages the 1/Rp
contribution to Re Z), and (a, n) from a power-law refit to the mean of the four
a w^n curves over 1-100 MHz, since the individual (a, n) pairs (a spans 26 orders
of magnitude) cannot be averaged term by term. At <= 1 MHz all of these are
negligible corrections to 1/(w C).

## Final parameter table

Contents of `fits/preamp_fits.py` (C in pF, L in nH, Rs in Ohm, a in Ohm s^n, Rp in kOhm).

Spare preamps (BNL flatsat, Bode 500, from Toki's presentation):

| preamp | C (pF) | L (nH) | Rs (Ohm) | n    | a (Ohm s^n) | Rp (kOhm) |
|--------|--------|--------|----------|------|-------------|-----------|
| fmpre0 | 33.79  | 35.0   | 0.21     | 1.39 | 5.2e-12     | 1101      |
| fmpre2 | 35.32  | 36.4   | 2.16     | 2.26 | 6.9e-20     | 331       |
| fmpre5 | 34.89  | 34.7   | 3.71     | 4.35 | 2.2e-38     | 288       |
| fmpre7 | 36.85  | 38.5   | 1.65     | 1.88 | 1.9e-16     | 375       |

Flight preamps (SDL, LCX200; C from the 1 MHz point, rest = spare averages):

| channel | preamp | C (pF)        | L (nH) | Rs (Ohm) | n    | a (Ohm s^n) | Rp (kOhm) |
|---------|--------|---------------|--------|----------|------|-------------|-----------|
| CH0     | fmpre6 | 33.12 +- 0.02 | 36.2   | 1.93     | 1.56 | 8.1e-14     | 397       |
| CH1     | fmpre3 | 33.29 +- 0.11 | 36.2   | 1.93     | 1.56 | 8.1e-14     | 397       |
| CH2     | fmpre4 | 33.44 +- 0.07 | 36.2   | 1.93     | 1.56 | 8.1e-14     | 397       |
| CH3     | fmpre1 | 33.17 +- 0.06 | 36.2   | 1.93     | 1.56 | 8.1e-14     | 397       |

The C errors are statistical only (standard error of the 5 repeats at 1 MHz);
allow a systematic of about +- 2 pF, see above.

To regenerate everything:

    python3 fits/fit_bode_spares.py
    python3 fits/fit_lcx200_flight.py
