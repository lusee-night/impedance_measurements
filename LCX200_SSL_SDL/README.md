# Preamp input impedance: LCX200 LCR meter at SDL/SSL

Summary of the email discussion (Feb–Mar 2026). Source threads:
"pre-amp impedance measurements full LuSEE freq range" (23 Feb – 9 Mar 2026)
and "Analyzing JFET impedance measurements" (25–27 Mar 2026).

## People
Kaja Rotermund (LBL, analysis), Emi Tamura, Jillian Donohue, Kayla Hernandez,
Sven Herrmann, Paul O'Connor (BNL), Ryan McLean (SSL), Toki Suzuki (LBL),
Stuart Bale (SSL/UCB), Anže Slosar.

## What was measured
- Preamp input impedance measured with an LCX200 LCR meter at SDL, by Jillian
  Donohue and Ryan McLean with the SSL team.
- Frequencies: 100, 180, 320, 560, 1000 kHz. Limited to 1 MHz by the
  instrument.
- Units in the spreadsheet: frequency in kHz, impedance in kOhm, phase in
  degrees (confirmed by Emi, 25 Feb).
- Data were written down by hand. The original Excel sheet could not be found
  (lost to Slack paywalling); Emi later dug up a copy.
- Measurement procedure: LST-SPT-0110, "Measurement_preamp-input-impedance_rev0"
  (PDF attached to Jillian's email of 9 Mar 2026).
- Per Kaja, 5 measurements were made per preamp at each frequency.
- It is not stated in the emails whether these were the flight-model preamps
  themselves.

## Gaps
- Kaja wanted the full 100 kHz – 50 MHz range to set uncertainties on the
  voltage power spectrum. VNA measurements of 1–50 MHz were expected but do
  not seem to exist. Jillian: nothing above 1 MHz was taken.
- This gap motivated the later Bode 500 measurements at BNL
  (see `../BODE500_BNL/`).

## Analysis
- Kaja (25 Mar): fitted a linear regression to resistance vs frequency in
  log-log space; the derived input capacitance was frequency dependent. Asked
  how to treat that.
- Emi (27 Mar): Kaja's fit was not acceptable, since data errors were not
  considered.
  - Series R + C model, Z = R + 1/(jwC): real part fine, imaginary part not
    statistically acceptable.
  - Adding series inductance, Z = R + jwL + 1/(jwC), or a parallel resistance
    gives reasonable fits, with C = 37 pF and 40 pF respectively.
  - Unclear which model is physical, so Emi recommends input capacitance
    **39 +- 2 pF for ch0**. Only ch0 was shown.
