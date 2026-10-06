# Data sources

Both files in this directory are **real** data, shipped verbatim from the GTEx
Consortium's v8 annotation bucket. No value in either file was simulated,
resampled, filtered or edited; only the filenames changed.

| file | upstream | rows x columns | sha256 |
|---|---|---|---|
| `gtex_v8_sample_attributes.tsv` | `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` | 22,951 x 63 | `74f6ab4c34ed2648d708a0ae6e6dff324f6c86ea723ae7d1c37d76f5221148f0` |
| `gtex_v8_subject_phenotypes.tsv` | `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | 980 x 4 | `821bdaff39e7a9a1d166919b3c786724c2b79c2861aeb936a2537a0f59b066f7` |

Both were downloaded from
`https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/`.

- `gtex_v8_sample_attributes.tsv` is one row per assayed GTEx sample. `SAMPID` is
  the sample key and encodes the donor as its first two `-`-separated fields;
  `SMTSD` is the sampled tissue site, `SMTS` its organ system, `SMAFRZE` the
  assay the sample was used for, `SMTSISCH` the post-mortem ischemic time in
  minutes (blank for some samples, negative for pre-mortem blood draws) and
  `SMRIN` the RNA integrity number. `SUBJID` is not a column in this file.
- `gtex_v8_subject_phenotypes.tsv` is one row per donor, keyed by `SUBJID`, with
  `SEX`, `AGE` as a ten-year bracket string (for example `60-69`) and `DTHHRDY`,
  the 0-4 Hardy death classification (blank for 19 donors).

## Provenance

- Donor: GTEx Consortium, *The GTEx Consortium atlas of genetic regulatory
  effects across human tissues*, Science 2020, DOI `10.1126/science.aaz1776`.
- Catalogue entry: `gtex-v8-sample-attributes` in
  `docs/sources/verified-quantitative.json` (its url is the sample-attributes
  file; its note names the subject-phenotype companion on the same bucket).
- Licence: GTEx open access. Summary-level annotation files such as these are
  downloadable without registration; raw sequence data require dbGaP
  authorisation and are not used here.

The machine-readable record of the same two files, with digests, is
`authoring/provenance/sources.json`; `authoring/provenance/fetch.py` re-downloads
them and refuses to rebuild the fixture if the upstream bytes have moved.
