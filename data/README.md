# Data

All data used in the paper are open. Nothing here needs credentialed access.

| Folder | Source | Licence | What is included |
| --- | --- | --- | --- |
| `eicu-demo/` | eICU Collaborative Research Database Demo v2.0.1 (PhysioNet) | Open Data Commons Open Database License v1.0 (ODbL) | The six tables used, unchanged, with the original licence and checksums |
| `nhanes/` | National Health and Nutrition Examination Survey, 2011–2018 (CDC/NCHS) | Public-use files | A download script and checksums of the files used |
| `derived/` | Analysis files written by `paper/analysis/eicu_aki.R` and `paper/analysis/nhanes.R` | eICU file: ODbL; NHANES file: public use | One row per patient: model risk, recorded outcome, reference outcome |

## eICU demo

Johnson A, Pollard T, Badawi O, Raffa J. eICU Collaborative Research Database Demo (version 2.0.1). PhysioNet; 2021. <https://doi.org/10.13026/4mxk-na84>

Please also cite PhysioNet as its website asks. The tables in `eicu-demo/` are redistributed under the ODbL; any database made from them must also be shared under the ODbL.

## NHANES

Files are downloaded from <https://wwwn.cdc.gov/nchs/nhanes/> by running, from the repository root,

```sh
Rscript data/nhanes/download_nhanes.R
```

`SHA256SUMS.txt` lists the exact files used. Users of NHANES data agree to the NCHS data use restrictions, including never trying to identify any participant.
