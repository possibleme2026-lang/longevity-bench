# Data sources

**Everything in this directory is REAL measured data.** Every number in
`nsc740_viability.csv` is a growth-inhibition endpoint that the NIH National
Cancer Institute's Developmental Therapeutics Program (NCI DTP) measured for
compound **NSC 740 -- methotrexate** in a human tumour cell line. Nothing here is
simulated, and no value has been imputed, rescaled or rounded.

## `nsc740_viability.csv`

One row per (cell line, concentration): 300 measurements = 60 cell lines x 5
concentrations. Columns
`cell_line,panel,concentration_log10m,endpoint_value,dispersion,n_experiments`.

| column | meaning |
|---|---|
| `cell_line` | the NCI-60 line the measurement belongs to |
| `panel` | the tissue panel DTP assigns that line to |
| `concentration_log10m` | compound concentration, log10 mol/L |
| `endpoint_value` | DTP's growth endpoint at that concentration: 100 is growth equal to the untreated control, 0 is no net growth over the assay, negative values are net cell loss |
| `dispersion` | DTP's spread of that endpoint across the replicate experiments behind it |
| `n_experiments` | how many replicate experiments the value averages (3 or 4) |

It is derived from the NCI DTP dose-response servlet for NSC 740:

```
https://dtp.cancer.gov/dtpstandard/servlet/SubSets?subsettype=DOSE&searchtype=NSC&searchlist=740&outputformat=xml
```

a 207,382-byte XML document (`sha256
207c650f9bea5ebb30ca1d4d3dfa0c3dd91a10520f8ce540b759e7a4baae543f`), public
domain, no registration. Catalogue label `nci-dtp-dose-response-nsc740`; the
catalogue entry's PMID 16990858 is the NCI-60 screen publication (Shoemaker,
2006), not the source of these numbers -- the servlet is. The URL, the byte count
and that SHA256 are recorded in `authoring/provenance/fetch.py`, which
re-downloads the payload, parses it and refuses to continue if either the
upstream digest or the rebuilt fixture digest has moved. The committed
derivation is `authoring/provenance/build_fixture.py`. Neither script is mounted
into the agent image or the verifier image.

### What the derivation did

The upstream document carries 133 `<CellInhibitionDoseResponse>` records, one per
cell line, and 1,724 `<DoseDataPoint>` measurements in total. Each record merges
every concentration series DTP ran for that line on this compound, and those
series come from different screening campaigns that are not comparable:

| campaign | concentrations (log10 M) | experiments per point | lines reporting it |
|---|---|---|---|
| A | `-9.0, -8.0, -7.0, -6.0, -5.0` (ten-fold) | 3 or 4 | 60 |
| B | `-7.6, -6.6, -5.6, -4.6, -3.6` (ten-fold) | 24 to 27 | 125 |
| C | `-6.31, -6.01, -5.71, -5.41, -5.11, -4.81, -4.51, -4.2, -3.9` (about two-fold) | 2 | 63 |
| -- | single extra concentrations | 1, `Dispersion = 0.0` | 89 records |

`build_fixture.py` keeps campaign **A** and nothing else:

1. keep the 60 cell lines whose record reports **all five** concentrations of
   campaign A;
2. keep those five measurements per line and drop the record's campaign-B,
   campaign-C and single-experiment points;
3. write `cell_line`, `panel`, `concentration_log10m`, `endpoint_value`,
   `dispersion` and `n_experiments` as they appear upstream.

No kept value is altered. Campaign A is the one grid that is simultaneously
(a) complete on a large contiguous block of cell lines, (b) uniform, so the
concentration axis is not a mixture of ten-fold and two-fold spacing, and
(c) replicated three or four times at every point, so every measurement has a
usable standard error -- 89 records carry single-experiment points whose
reported dispersion is `0.0`, which cannot be weighted at all. Mixing campaigns
would also let replicate weight, not chemistry, drive the answer: campaign B's
points average 24-27 experiments against campaign A's 3-4.

Campaign A is also the only campaign that brackets this compound's potency. On
campaign B's grid, 32 of the 60 lines -- including 32 of the 45 whose GI50
campaign A resolves -- are already inhibited past 50 at the lowest concentration
tested, so they would be censored rather than measured; campaign A leaves 15 of
60 censored instead of 41.

### What is in it

* 60 cell lines over 9 panels: Colon (7), Melanoma (10), Non-Small Cell Lung
  (9), Renal (8), Ovarian (7), Central Nervous System (6), Leukemia (6), Breast
  (5) and Prostate (2).
* 5 concentrations per line, 1 nM to 10 uM in ten-fold steps.
* `n_experiments` is 4 for 228 rows and 3 for 72 rows; `dispersion` is strictly
  positive in every row, so every measurement has a usable standard error.
* 45 of the 60 lines take the endpoint to 50 or below somewhere in the tested
  range; for 15 lines the endpoint never falls to 50, so their GI50 lies above
  10 uM and is not measurable in this screen. That asymmetry is a property of
  methotrexate in the NCI-60, not of the derivation: the sensitive lines are the
  leukemias (K-562, SR, CCRF-CEM), and the resistant ones are spread across
  several panels.

### What is deliberately not here

The file holds no derived quantity. There is no GI50, no standard error, no
panel mean and no ranking -- those are what the task asks for.
