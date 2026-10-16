# Preamp input impedance: Bode 500 at BNL (flatsat)

Summary of the email discussion (Apr–Jun 2026). Source threads: "Send BODE
system to BNL" (15–20 Apr) and "flat sat" (19 May – 25 Jun).

## People
Anže Slosar, Emi Tamura (measurements), Toki Suzuki (LBL, owner of the Bode
and analysis), Paul O'Connor, Sven Herrmann.

## Setup
- Bode 500 (Omicron) impedance analyzer, shipped from LBL to BNL in April.
  Manual: Omicron Bode Analyzer Suite user manual.
- Goal: input impedance of all 4 preamps across the full band (up to 50 MHz).
- Calibration: open/short/load at the preamp box input connector, to de-embed
  the cable (Paul's advice). One-port mode.
- Anže's own preliminary measurement suggested the input stops being purely
  capacitive above about 40 MHz.
- Flatsat: a single flatsat in the cleanroom with 4 preamps. CH2 (counting
  from 0) is on the AWG, the others are terminated with 50 ohm. Working,
  including AWG, as of 25 Jun.

## Results
- **FM preamp #0 (Emi, 28 May):** purely capacitive up to 10 MHz. Best-fit
  C = 33.8 pF, consistent with an earlier 34.9 pF. Measured up to 100 MHz;
  above 10 MHz the impedance is lower than the capacitive model, so another
  component is needed. Averaged data of 10 acquisitions are in
  `fmpre0_impedance.dat` (frequency, average, sd). Emi has the raw complex
  data.
- **Cross-talk (Emi, 3–4 Jun):** S12 and S21 between #0 (port 1) and #2
  (port 2), stimulus -25 dBm, with shared and with separate power supplies.
  S21 is lower than S12 in both cases. Anže reads the off-diagonal terms of
  Z_L as at least 50 dB below the diagonal, hence negligible. Toki: coupling
  is present but at a low level.
- Caveat: the flatsat has unboxed or 3D-printed enclosures, different
  cabling, and no carousel, so it is not representative of the flight
  electrical configuration. Toki would like four FM preamps mounted on a
  metal plate. Emi recommends measuring the actual FMs on LuSEE-Night.
- **Toki's analysis (9 Jun):** sent to the thread; Toki agreed to present it at
  the antenna meeting. Includes a 35 nH series inductance.
- **Paul's open questions (9 Jun):**
  1. Does the 35 nH series inductance improve or worsen the high-frequency
     antenna/preamp mismatch, or add resonances?
  2. Gate resistance will shift at the temperature extremes. How important is
     the temperature dependence of Zin?
  Answers were not seen.

## Status
- Bode was to be returned to LBNL after the measurements (Emi asked on 4 Jun;
  Toki gave a shipping address).
- Not yet measured as of these emails: the four EM preamps in the cleanroom
  (Emi offered).

## Related
- Earlier LCR-meter measurements (to 1 MHz): `../LCX200_SSL_SDL/`.
