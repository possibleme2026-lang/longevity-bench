# Data sources

**Everything in this directory is REAL measured data.** Every number in
`itp_lifespan.csv` is an observation recorded by the NIH Interventions Testing
Program (ITP) on a UM-HET3 mouse, downloaded from the Mouse Phenome Database
(MPD) at The Jackson Laboratory. Nothing here is simulated, and no value has
been imputed, rescaled or rounded.

## `itp_lifespan.csv`

One row per animal: 5,255 mice. Columns
`cohort,mouse_id,site,sex,group,dose_ppm,age_initiation_months,dead,age_days`.

It is the union of two MPD curation files, both free public downloads under the
Mouse Phenome Database Terms of Use (The Jackson Laboratory):

| MPD curation id | workbook | catalogue label | reference | mice |
|---|---|---|---|---|
| 772 | `ITP_C2021_Lifespan.xlsx` | `mpd-itp1-lifespan-cohort-2021-dose-levels` | PMID 19587680 | 2,388 |
| 707 | `ITP_C2020_Lifespan.xlsx` | `mpd-itp1-lifespan-cohort-2020-dose-levels` | PMID 27312235 | 2,869 |

Both URLs are recorded in `authoring/provenance/fetch.py`, which re-downloads
them, checks the ZIP magic bytes and the `xl/workbook.xml` member (an HTTP 200
on its own is not evidence that a data file arrived), refuses to continue if
either payload's SHA256 has moved, and rebuilds this file byte for byte. The
committed derivation is `authoring/provenance/build_fixture.py`. Neither script
is mounted into the agent image or the verifier image.

The two PMIDs are the ITP programme papers the catalogue cites for these
cohorts; PMID 19587680 (Harrison et al. 2009, rapamycin) is the founding ITP
lifespan report and PMID 27312235 (Strong et al. 2016) is a later multi-arm
report. Neither paper is the source of these numbers -- the MPD curation files
are. The drug codes below are the ITP's own arm labels.

### What the derivation did

Only these changes, in this order:

1. kept every source column except the constant `population` (= `UM-HET3`) and
   `status` (which duplicates `dead` as text);
2. renamed the columns to one schema across the two workbooks, whose headers
   genuinely differ (`dose` vs `Rx(ppm)`, `age` vs `age(days)`,
   `age_initiation` vs `age_initiation(mo)`);
3. turned the 2020 workbook's literal `NA` age at initiation into an empty
   field, matching the 2021 workbook, which leaves it blank;
4. parsed `dose` to a number;
5. dropped the two byte-identical repeat records of animal `JL21949` that the
   2020 workbook carries three times.

No animal is excluded on the basis of its age, weight or survival time. The
shortest recorded lifetime is 22 days and the longest 1,461 days; both are real
records and both are kept.

### What is in it

* Two cohorts, `C2021` (2,388 mice) and `C2020` (2,867 mice), tested at three
  sites (`TJL`, `UM`, `UT`).
* Both sexes: 2,488 females and 2,767 males.
* 16 distinct `group` labels: `Control` plus 15 dosed arms. `Control` animals
  carry `dose_ppm` `0` and an empty `age_initiation_months` because they were
  never dosed -- all 1,169 blank initiation ages are `Control` animals, and no
  dosed animal has one.
* 9 distinct non-zero dose levels -- 0.6, 3, 5, 30, 60, 180, 300, 10,000 and
  20,000 ppm -- over 34 cohort-by-sex-by-arm strata and 28 cohort-by-sex-by-dose
  dose levels. Three dose levels pool more than one arm: 30 ppm in C2021 pools
  2BA and DCA, two different compounds; 30 ppm and 180 ppm in C2020 each pool
  the two initiation ages of one compound (HYD_6 with HYD_16, and Cana_6 with
  Cana_16).
* `dead` = 1 for 5,122 animals and 0 for 133: those 133 were removed from the
  study alive and are right-censored at `age_days`.

### What is deliberately not here

The file holds no derived quantity. There is no lifespan ratio, no p-value, no
survival estimate, and no ranking -- those are what the task asks for.
