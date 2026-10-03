# Table S-dsep: d-separation test of both DAG variants against data

Task 4.3. One row per `impliedConditionalIndependencies(dag, type="missing.edge")`
test, for both DAG variants from 03_dag.R. The compositional `Benthos`/
`Benthos_lag` nodes are tested as their full 3-dimensional ilr vector via
`cis.pillai` (canonical correlation), not reduced to a scalar. Tests
involving Rd/Pmax/DayTemp/Flow/Season use `panel_lter1` (LTER_1 only,
n~17-18, PROXY Rd/Pmax from pp_day.csv -- NOT the Section 6 PI-curve
model); all other tests use `panel_full` (6 sites, n~107-120).
`p_adj_BH` is the Benjamini-Hochberg-adjusted p-value within each DAG
variant (440 testable tests total across both variants). See
`qc_dsep.md` for the Shipley/Fisher's C global test and interpretation.

**Summary: 442 rows, 440 testable, 15 violated (BH-adjusted p<0.05), 2 not testable.**

| DAG variant | X | Y | Z | estimate | p | p (BH) | n | status |
|---|---|---|---|---|---|---|---|---|
| Rd -> N | Benthos | DayTemp | {DHW} | 0.416 | 0.43 | 0.78 | 18 | consistent with DAG |
| Rd -> N | Benthos | Flow | {} | 0.262 | 0.79 | 0.97 | 18 | consistent with DAG |
| Rd -> N | Benthos | N | {Corall, Herb, OtherFish, Rd, Time} | 0.703 | 0.027 | 0.26 | 17 | consistent with DAG |
| Rd -> N | Benthos | Season | {} | 0.322 | 0.66 | 0.95 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Corall | {Benthos} | 0.193 | 0.24 | 0.59 | 114 | consistent with DAG |
| Rd -> N | Benthos_lag | Cyclone | {} | 0.244 | 0.079 | 0.43 | 114 | consistent with DAG |
| Rd -> N | Benthos_lag | DHW | {Time} | 0.179 | 0.31 | 0.65 | 114 | consistent with DAG |
| Rd -> N | Benthos_lag | DayTemp | {DHW} | 0.504 | 0.24 | 0.59 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | DayTemp | {Time} | 0.477 | 0.29 | 0.63 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Flow | {} | 0.378 | 0.53 | 0.87 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Herb_lag | {Time} | 0.372 | 0.0009 | 0.034 | 114 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | Benthos_lag | N | {Corall, Herb, OtherFish, Rd, Time} | 0.583 | 0.13 | 0.52 | 17 | consistent with DAG |
| Rd -> N | Benthos_lag | N | {Benthos, DayTemp, Flow, Herb, Time} | 0.626 | 0.082 | 0.43 | 17 | consistent with DAG |
| Rd -> N | Benthos_lag | N | {Benthos, DayTemp, Herb, Season, Time} | 0.683 | 0.037 | 0.29 | 17 | consistent with DAG |
| Rd -> N | Benthos_lag | N | {Benthos, DHW, Herb, Time} | 0.103 | 0.77 | 0.97 | 107 | consistent with DAG |
| Rd -> N | Benthos_lag | N_lag | {Time} | 0.110 | 0.74 | 0.96 | 107 | consistent with DAG |
| Rd -> N | Benthos_lag | OtherFish | {Benthos, Time} | 0.161 | 0.4 | 0.76 | 114 | consistent with DAG |
| Rd -> N | Benthos_lag | Pmax | {Benthos, DayTemp, Flow, N} | 0.687 | 0.035 | 0.29 | 17 | consistent with DAG |
| Rd -> N | Benthos_lag | Pmax | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | NA | NA | NA | 17 | not testable |
| Rd -> N | Benthos_lag | Pmax | {Benthos, Corall, DHW, Herb, N, OtherFish, Rd} | NA | NA | NA | 17 | not testable |
| Rd -> N | Benthos_lag | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.719 | 0.015 | 0.19 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Pmax | {Benthos, DayTemp, Herb, Season, Time} | 0.727 | 0.013 | 0.19 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Pmax | {Benthos, DHW, Herb, Time} | 0.574 | 0.12 | 0.51 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.655 | 0.044 | 0.32 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.665 | 0.038 | 0.29 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Rd | {Benthos, DHW, Herb, OtherFish} | 0.686 | 0.027 | 0.26 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.723 | 0.013 | 0.19 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.732 | 0.011 | 0.19 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Rd | {Benthos, DHW, Herb, Time} | 0.738 | 0.0099 | 0.19 | 18 | consistent with DAG |
| Rd -> N | Benthos_lag | Season | {} | 0.357 | 0.58 | 0.9 | 18 | consistent with DAG |
| Rd -> N | COTS | Corall | {Benthos} | 0.157 | 0.087 | 0.43 | 120 | consistent with DAG |
| Rd -> N | COTS | Cyclone | {} | 0.002 | 0.98 | 1 | 120 | consistent with DAG |
| Rd -> N | COTS | DHW | {Time} | -0.091 | 0.32 | 0.66 | 120 | consistent with DAG |
| Rd -> N | COTS | DayTemp | {DHW} | 0.287 | 0.25 | 0.59 | 18 | consistent with DAG |
| Rd -> N | COTS | DayTemp | {Time} | 0.322 | 0.19 | 0.57 | 18 | consistent with DAG |
| Rd -> N | COTS | Flow | {} | 0.319 | 0.2 | 0.57 | 18 | consistent with DAG |
| Rd -> N | COTS | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | 0.093 | 0.33 | 0.66 | 114 | consistent with DAG |
| Rd -> N | COTS | Herb_lag | {Time} | -0.102 | 0.28 | 0.63 | 114 | consistent with DAG |
| Rd -> N | COTS | N | {Corall, Herb, OtherFish, Rd, Time} | -0.090 | 0.73 | 0.96 | 17 | consistent with DAG |
| Rd -> N | COTS | N | {Benthos, DayTemp, Flow, Herb, Time} | 0.008 | 0.98 | 1 | 17 | consistent with DAG |
| Rd -> N | COTS | N | {Benthos, DayTemp, Herb, Season, Time} | -0.125 | 0.63 | 0.94 | 17 | consistent with DAG |
| Rd -> N | COTS | N | {Benthos, DHW, Herb, Time} | 0.322 | 0.00072 | 0.032 | 107 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | COTS | N | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.445 | 0.074 | 0.43 | 17 | consistent with DAG |
| Rd -> N | COTS | N | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.334 | 0.19 | 0.57 | 17 | consistent with DAG |
| Rd -> N | COTS | N | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.370 | 8.8e-05 | 0.0068 | 107 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | COTS | N_lag | {Time} | 0.495 | 5.9e-08 | 1.3e-05 | 107 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | COTS | OtherFish | {Benthos, Time} | 0.018 | 0.85 | 0.97 | 120 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, DayTemp, Flow, N} | -0.411 | 0.1 | 0.46 | 17 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | -0.705 | 0.0016 | 0.045 | 17 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | COTS | Pmax | {Benthos, Corall, DHW, Herb, N, OtherFish, Rd} | -0.807 | 9e-05 | 0.0068 | 17 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | COTS | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | -0.456 | 0.057 | 0.36 | 18 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, DayTemp, Herb, Season, Time} | -0.369 | 0.13 | 0.52 | 18 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, DHW, Herb, Time} | -0.568 | 0.014 | 0.19 | 18 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.484 | 0.042 | 0.31 | 18 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.470 | 0.049 | 0.34 | 18 | consistent with DAG |
| Rd -> N | COTS | Pmax | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.703 | 0.0011 | 0.037 | 18 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | COTS | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.137 | 0.59 | 0.9 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.079 | 0.75 | 0.96 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, DHW, Herb, OtherFish} | 0.197 | 0.43 | 0.78 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.054 | 0.83 | 0.97 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.041 | 0.87 | 0.97 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, DHW, Herb, Time} | 0.188 | 0.45 | 0.79 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | 0.038 | 0.88 | 0.97 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | 0.050 | 0.84 | 0.97 | 18 | consistent with DAG |
| Rd -> N | COTS | Rd | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.168 | 0.51 | 0.86 | 18 | consistent with DAG |
| Rd -> N | COTS | Season | {} | 0.367 | 0.13 | 0.52 | 18 | consistent with DAG |
| Rd -> N | Corall | Cyclone | {Benthos} | -0.093 | 0.31 | 0.65 | 120 | consistent with DAG |
| Rd -> N | Corall | DHW | {Benthos} | 0.000 | 1 | 1 | 120 | consistent with DAG |
| Rd -> N | Corall | DayTemp | {DHW} | 0.049 | 0.85 | 0.97 | 18 | consistent with DAG |
| Rd -> N | Corall | DayTemp | {Benthos} | 0.348 | 0.16 | 0.54 | 18 | consistent with DAG |
| Rd -> N | Corall | Flow | {} | 0.021 | 0.94 | 1 | 18 | consistent with DAG |
| Rd -> N | Corall | Herb | {Benthos} | 0.113 | 0.22 | 0.58 | 120 | consistent with DAG |
| Rd -> N | Corall | Herb_lag | {Benthos} | -0.112 | 0.23 | 0.59 | 114 | consistent with DAG |
| Rd -> N | Corall | N_lag | {Benthos} | 0.022 | 0.82 | 0.97 | 107 | consistent with DAG |
| Rd -> N | Corall | OtherFish | {Benthos} | 0.157 | 0.087 | 0.43 | 120 | consistent with DAG |
| Rd -> N | Corall | Pmax | {Benthos, DayTemp, Flow, N} | 0.125 | 0.63 | 0.94 | 17 | consistent with DAG |
| Rd -> N | Corall | Season | {} | -0.118 | 0.64 | 0.94 | 18 | consistent with DAG |
| Rd -> N | Corall | Time | {Benthos} | -0.134 | 0.14 | 0.53 | 120 | consistent with DAG |
| Rd -> N | Cyclone | DHW | {} | -0.110 | 0.23 | 0.59 | 120 | consistent with DAG |
| Rd -> N | Cyclone | DayTemp | {} | -0.630 | 0.005 | 0.11 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Flow | {} | -0.003 | 0.99 | 1 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | -0.073 | 0.44 | 0.78 | 114 | consistent with DAG |
| Rd -> N | Cyclone | Herb_lag | {} | -0.179 | 0.057 | 0.36 | 114 | consistent with DAG |
| Rd -> N | Cyclone | N | {Corall, Herb, OtherFish, Rd, Time} | -0.250 | 0.33 | 0.66 | 17 | consistent with DAG |
| Rd -> N | Cyclone | N | {Benthos, DayTemp, Flow, Herb, Time} | -0.317 | 0.21 | 0.58 | 17 | consistent with DAG |
| Rd -> N | Cyclone | N | {Benthos, DayTemp, Herb, Season, Time} | -0.400 | 0.11 | 0.49 | 17 | consistent with DAG |
| Rd -> N | Cyclone | N | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.105 | 0.69 | 0.96 | 17 | consistent with DAG |
| Rd -> N | Cyclone | N | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.086 | 0.74 | 0.96 | 17 | consistent with DAG |
| Rd -> N | Cyclone | N | {Benthos, DHW, Herb, Time} | 0.025 | 0.8 | 0.97 | 107 | consistent with DAG |
| Rd -> N | Cyclone | N | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.054 | 0.58 | 0.9 | 107 | consistent with DAG |
| Rd -> N | Cyclone | N_lag | {} | 0.051 | 0.6 | 0.92 | 107 | consistent with DAG |
| Rd -> N | Cyclone | OtherFish | {Benthos, Time} | -0.210 | 0.021 | 0.24 | 120 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, DayTemp, Flow, N} | -0.053 | 0.84 | 0.97 | 17 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | -0.140 | 0.59 | 0.9 | 17 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.004 | 0.99 | 1 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, DayTemp, Herb, Season, Time} | -0.019 | 0.94 | 1 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.053 | 0.83 | 0.97 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.140 | 0.58 | 0.9 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, Corall, DHW, Herb, N, OtherFish, Rd} | 0.568 | 0.017 | 0.21 | 17 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, DHW, Herb, Time} | 0.384 | 0.12 | 0.5 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Pmax | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.402 | 0.098 | 0.45 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | -0.354 | 0.15 | 0.53 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | -0.379 | 0.12 | 0.51 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, DayTemp, Flow, Herb, Time} | -0.271 | 0.28 | 0.63 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, DayTemp, Herb, Season, Time} | -0.267 | 0.28 | 0.63 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.406 | 0.095 | 0.45 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.429 | 0.075 | 0.43 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, DHW, Herb, OtherFish} | -0.432 | 0.074 | 0.43 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, DHW, Herb, Time} | -0.356 | 0.15 | 0.53 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Rd | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.554 | 0.017 | 0.21 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Season | {} | -0.367 | 0.13 | 0.52 | 18 | consistent with DAG |
| Rd -> N | Cyclone | Time | {} | -0.179 | 0.05 | 0.34 | 120 | consistent with DAG |
| Rd -> N | DHW | Flow | {} | 0.045 | 0.86 | 0.97 | 18 | consistent with DAG |
| Rd -> N | DHW | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | -0.006 | 0.95 | 1 | 114 | consistent with DAG |
| Rd -> N | DHW | Herb_lag | {Time} | -0.014 | 0.89 | 0.98 | 114 | consistent with DAG |
| Rd -> N | DHW | N | {Corall, Herb, OtherFish, Rd, Time} | -0.010 | 0.97 | 1 | 17 | consistent with DAG |
| Rd -> N | DHW | N | {Benthos, DayTemp, Flow, Herb, Time} | 0.210 | 0.42 | 0.77 | 17 | consistent with DAG |
| Rd -> N | DHW | N | {Benthos, DayTemp, Herb, Season, Time} | 0.308 | 0.23 | 0.59 | 17 | consistent with DAG |
| Rd -> N | DHW | N | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.482 | 0.05 | 0.34 | 17 | consistent with DAG |
| Rd -> N | DHW | N | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.516 | 0.034 | 0.29 | 17 | consistent with DAG |
| Rd -> N | DHW | N_lag | {Time} | -0.026 | 0.79 | 0.97 | 107 | consistent with DAG |
| Rd -> N | DHW | OtherFish | {Benthos, Time} | 0.031 | 0.74 | 0.96 | 120 | consistent with DAG |
| Rd -> N | DHW | Pmax | {Benthos, DayTemp, Flow, N} | 0.134 | 0.61 | 0.92 | 17 | consistent with DAG |
| Rd -> N | DHW | Pmax | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | -0.000 | 1 | 1 | 17 | consistent with DAG |
| Rd -> N | DHW | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.056 | 0.83 | 0.97 | 18 | consistent with DAG |
| Rd -> N | DHW | Pmax | {Benthos, DayTemp, Herb, Season, Time} | 0.171 | 0.5 | 0.85 | 18 | consistent with DAG |
| Rd -> N | DHW | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.332 | 0.18 | 0.57 | 18 | consistent with DAG |
| Rd -> N | DHW | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.295 | 0.23 | 0.59 | 18 | consistent with DAG |
| Rd -> N | DHW | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.242 | 0.33 | 0.66 | 18 | consistent with DAG |
| Rd -> N | DHW | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.314 | 0.2 | 0.57 | 18 | consistent with DAG |
| Rd -> N | DHW | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.269 | 0.28 | 0.63 | 18 | consistent with DAG |
| Rd -> N | DHW | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.321 | 0.19 | 0.57 | 18 | consistent with DAG |
| Rd -> N | DHW | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.285 | 0.25 | 0.59 | 18 | consistent with DAG |
| Rd -> N | DHW | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.306 | 0.22 | 0.58 | 18 | consistent with DAG |
| Rd -> N | DHW | Season | {} | 0.502 | 0.034 | 0.29 | 18 | consistent with DAG |
| Rd -> N | DayTemp | Flow | {Season} | 0.097 | 0.7 | 0.96 | 18 | consistent with DAG |
| Rd -> N | DayTemp | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | 0.202 | 0.42 | 0.77 | 18 | consistent with DAG |
| Rd -> N | DayTemp | Herb | {DHW} | 0.282 | 0.26 | 0.6 | 18 | consistent with DAG |
| Rd -> N | DayTemp | Herb_lag | {Time} | 0.423 | 0.08 | 0.43 | 18 | consistent with DAG |
| Rd -> N | DayTemp | Herb_lag | {DHW} | 0.359 | 0.14 | 0.53 | 18 | consistent with DAG |
| Rd -> N | DayTemp | N | {Corall, Herb, OtherFish, Rd, Time} | -0.158 | 0.55 | 0.89 | 17 | consistent with DAG |
| Rd -> N | DayTemp | N | {Benthos, Corall, DHW, Herb, OtherFish, Rd} | -0.112 | 0.67 | 0.95 | 17 | consistent with DAG |
| Rd -> N | DayTemp | N_lag | {Time} | 0.358 | 0.14 | 0.53 | 18 | consistent with DAG |
| Rd -> N | DayTemp | N_lag | {DHW} | 0.030 | 0.9 | 0.98 | 18 | consistent with DAG |
| Rd -> N | DayTemp | OtherFish | {Benthos, Time} | 0.240 | 0.34 | 0.66 | 18 | consistent with DAG |
| Rd -> N | DayTemp | OtherFish | {DHW} | 0.254 | 0.31 | 0.65 | 18 | consistent with DAG |
| Rd -> N | DayTemp | Time | {DHW} | 0.137 | 0.59 | 0.9 | 18 | consistent with DAG |
| Rd -> N | Flow | Herb | {} | 0.213 | 0.4 | 0.75 | 18 | consistent with DAG |
| Rd -> N | Flow | Herb_lag | {} | 0.274 | 0.27 | 0.62 | 18 | consistent with DAG |
| Rd -> N | Flow | N | {Corall, Herb, OtherFish, Rd, Time} | -0.239 | 0.36 | 0.68 | 17 | consistent with DAG |
| Rd -> N | Flow | N | {Benthos, Corall, DHW, Herb, OtherFish, Rd} | -0.045 | 0.86 | 0.97 | 17 | consistent with DAG |
| Rd -> N | Flow | N | {Benthos, Corall, DayTemp, Herb, OtherFish, Rd, Season} | -0.004 | 0.99 | 1 | 17 | consistent with DAG |
| Rd -> N | Flow | N_lag | {} | -0.073 | 0.77 | 0.97 | 18 | consistent with DAG |
| Rd -> N | Flow | OtherFish | {} | -0.230 | 0.36 | 0.68 | 18 | consistent with DAG |
| Rd -> N | Flow | Time | {} | 0.105 | 0.68 | 0.96 | 18 | consistent with DAG |
| Rd -> N | Herb | N_lag | {Benthos, Benthos_lag, Herb_lag, Time} | -0.114 | 0.24 | 0.59 | 107 | consistent with DAG |
| Rd -> N | Herb | OtherFish | {Benthos, Time} | 0.249 | 0.006 | 0.12 | 120 | consistent with DAG |
| Rd -> N | Herb | Pmax | {Benthos, DayTemp, Flow, N} | 0.066 | 0.8 | 0.97 | 17 | consistent with DAG |
| Rd -> N | Herb | Season | {} | 0.410 | 0.091 | 0.44 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | N | {Corall, Herb, OtherFish, Rd, Time} | -0.352 | 0.17 | 0.55 | 17 | consistent with DAG |
| Rd -> N | Herb_lag | N | {Benthos, DayTemp, Flow, Herb, Time} | 0.106 | 0.69 | 0.96 | 17 | consistent with DAG |
| Rd -> N | Herb_lag | N | {Benthos, DayTemp, Herb, Season, Time} | -0.024 | 0.93 | 1 | 17 | consistent with DAG |
| Rd -> N | Herb_lag | N | {Benthos, DHW, Herb, Time} | 0.273 | 0.0044 | 0.11 | 107 | consistent with DAG |
| Rd -> N | Herb_lag | N_lag | {Time} | 0.017 | 0.86 | 0.97 | 107 | consistent with DAG |
| Rd -> N | Herb_lag | OtherFish | {Benthos, Time} | -0.011 | 0.9 | 0.98 | 114 | consistent with DAG |
| Rd -> N | Herb_lag | Pmax | {Benthos, DayTemp, Flow, N} | 0.382 | 0.13 | 0.52 | 17 | consistent with DAG |
| Rd -> N | Herb_lag | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.292 | 0.24 | 0.59 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Pmax | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | -0.124 | 0.64 | 0.94 | 17 | consistent with DAG |
| Rd -> N | Herb_lag | Pmax | {Benthos, DayTemp, Herb, Season, Time} | 0.256 | 0.3 | 0.65 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Pmax | {Benthos, Corall, DHW, Herb, N, OtherFish, Rd} | -0.542 | 0.025 | 0.25 | 17 | consistent with DAG |
| Rd -> N | Herb_lag | Pmax | {Benthos, DHW, Herb, Time} | 0.159 | 0.53 | 0.87 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.228 | 0.36 | 0.69 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.287 | 0.25 | 0.59 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.240 | 0.34 | 0.66 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.312 | 0.21 | 0.57 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Rd | {Benthos, DHW, Herb, OtherFish} | 0.347 | 0.16 | 0.54 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Rd | {Benthos, DHW, Herb, Time} | 0.414 | 0.087 | 0.43 | 18 | consistent with DAG |
| Rd -> N | Herb_lag | Season | {} | 0.233 | 0.35 | 0.68 | 18 | consistent with DAG |
| Rd -> N | N | N_lag | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.046 | 0.65 | 0.95 | 101 | consistent with DAG |
| Rd -> N | N | N_lag | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.336 | 0.19 | 0.57 | 17 | consistent with DAG |
| Rd -> N | N | N_lag | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.203 | 0.44 | 0.78 | 17 | consistent with DAG |
| Rd -> N | N | N_lag | {Benthos, DHW, Herb, Time} | 0.023 | 0.82 | 0.97 | 101 | consistent with DAG |
| Rd -> N | N | N_lag | {Benthos, DayTemp, Herb, Season, Time} | 0.091 | 0.73 | 0.96 | 17 | consistent with DAG |
| Rd -> N | N | N_lag | {Benthos, DayTemp, Flow, Herb, Time} | 0.040 | 0.88 | 0.97 | 17 | consistent with DAG |
| Rd -> N | N | N_lag | {Corall, Herb, OtherFish, Rd, Time} | 0.117 | 0.65 | 0.95 | 17 | consistent with DAG |
| Rd -> N | N | Season | {DHW, DayTemp, Flow} | -0.045 | 0.86 | 0.97 | 17 | consistent with DAG |
| Rd -> N | N | Season | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.162 | 0.53 | 0.88 | 17 | consistent with DAG |
| Rd -> N | N | Season | {Benthos, DayTemp, Flow, Herb, Time} | 0.094 | 0.72 | 0.96 | 17 | consistent with DAG |
| Rd -> N | N | Season | {Benthos, Corall, DHW, Herb, OtherFish, Rd} | -0.207 | 0.43 | 0.78 | 17 | consistent with DAG |
| Rd -> N | N | Season | {Corall, Herb, OtherFish, Rd, Time} | -0.091 | 0.73 | 0.96 | 17 | consistent with DAG |
| Rd -> N | N_lag | OtherFish | {Benthos, Time} | 0.220 | 0.023 | 0.24 | 107 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, DayTemp, Flow, N} | -0.015 | 0.95 | 1 | 17 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.004 | 0.99 | 1 | 18 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | 0.093 | 0.71 | 0.96 | 18 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | 0.322 | 0.21 | 0.57 | 17 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, DayTemp, Herb, Season, Time} | -0.141 | 0.58 | 0.9 | 18 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | 0.007 | 0.98 | 1 | 18 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, Corall, DHW, Herb, N, OtherFish, Rd} | 0.018 | 0.95 | 1 | 17 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, DHW, Herb, Time} | -0.159 | 0.53 | 0.87 | 18 | consistent with DAG |
| Rd -> N | N_lag | Pmax | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.166 | 0.51 | 0.86 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | -0.070 | 0.78 | 0.97 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, DayTemp, Flow, Herb, Time} | -0.146 | 0.56 | 0.9 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.256 | 0.3 | 0.65 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | -0.101 | 0.69 | 0.96 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, DayTemp, Herb, Season, Time} | -0.175 | 0.49 | 0.84 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.341 | 0.17 | 0.55 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, DHW, Herb, OtherFish} | 0.094 | 0.71 | 0.96 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, DHW, Herb, Time} | 0.058 | 0.82 | 0.97 | 18 | consistent with DAG |
| Rd -> N | N_lag | Rd | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.341 | 0.17 | 0.55 | 18 | consistent with DAG |
| Rd -> N | N_lag | Season | {} | -0.324 | 0.19 | 0.57 | 18 | consistent with DAG |
| Rd -> N | OtherFish | Pmax | {Benthos, DayTemp, Flow, N} | -0.444 | 0.074 | 0.43 | 17 | consistent with DAG |
| Rd -> N | OtherFish | Season | {} | 0.096 | 0.7 | 0.96 | 18 | consistent with DAG |
| Rd -> N | Pmax | Rd | {Benthos, DayTemp, Flow, N} | 0.777 | 0.00024 | 0.014 | 17 | VIOLATED (BH-adjusted p < 0.05) |
| Rd -> N | Pmax | Season | {DHW, DayTemp, Flow} | -0.313 | 0.21 | 0.57 | 18 | consistent with DAG |
| Rd -> N | Pmax | Season | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.111 | 0.66 | 0.95 | 18 | consistent with DAG |
| Rd -> N | Pmax | Season | {Benthos, DayTemp, Flow, Herb, Time} | -0.187 | 0.46 | 0.79 | 18 | consistent with DAG |
| Rd -> N | Pmax | Season | {Benthos, DayTemp, Flow, N} | -0.146 | 0.58 | 0.9 | 17 | consistent with DAG |
| Rd -> N | Pmax | Time | {Benthos, Corall, DHW, Herb, N, OtherFish, Rd} | -0.273 | 0.29 | 0.63 | 17 | consistent with DAG |
| Rd -> N | Pmax | Time | {Benthos, Corall, DayTemp, Herb, N, OtherFish, Rd, Season} | -0.324 | 0.21 | 0.57 | 17 | consistent with DAG |
| Rd -> N | Pmax | Time | {Benthos, DayTemp, Flow, N} | 0.074 | 0.78 | 0.97 | 17 | consistent with DAG |
| Rd -> N | Rd | Season | {DHW, DayTemp, Flow} | -0.081 | 0.75 | 0.96 | 18 | consistent with DAG |
| Rd -> N | Rd | Season | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.088 | 0.73 | 0.96 | 18 | consistent with DAG |
| Rd -> N | Rd | Season | {Benthos, DayTemp, Flow, Herb, Time} | 0.021 | 0.93 | 1 | 18 | consistent with DAG |
| Rd -> N | Rd | Season | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.034 | 0.89 | 0.98 | 18 | consistent with DAG |
| Rd -> N | Rd | Time | {Benthos, DHW, Herb, OtherFish} | -0.042 | 0.87 | 0.97 | 18 | consistent with DAG |
| Rd -> N | Rd | Time | {Benthos, DayTemp, Herb, OtherFish, Season} | -0.025 | 0.92 | 1 | 18 | consistent with DAG |
| Rd -> N | Rd | Time | {Benthos, DayTemp, Flow, Herb, OtherFish} | -0.051 | 0.84 | 0.97 | 18 | consistent with DAG |
| Rd -> N | Season | Time | {} | 0.314 | 0.2 | 0.57 | 18 | consistent with DAG |
| Benthos -> N | Benthos | DayTemp | {DHW} | 0.416 | 0.43 | 0.78 | 18 | consistent with DAG |
| Benthos -> N | Benthos | Flow | {} | 0.262 | 0.79 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Benthos | Season | {} | 0.322 | 0.66 | 0.93 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Corall | {Benthos} | 0.193 | 0.24 | 0.59 | 114 | consistent with DAG |
| Benthos -> N | Benthos_lag | Cyclone | {} | 0.244 | 0.079 | 0.42 | 114 | consistent with DAG |
| Benthos -> N | Benthos_lag | DHW | {Time} | 0.179 | 0.31 | 0.65 | 114 | consistent with DAG |
| Benthos -> N | Benthos_lag | DayTemp | {DHW} | 0.504 | 0.24 | 0.59 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | DayTemp | {Time} | 0.477 | 0.29 | 0.64 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Flow | {} | 0.378 | 0.53 | 0.86 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Herb_lag | {Time} | 0.372 | 0.0009 | 0.039 | 114 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | Benthos_lag | N | {Benthos, Herb, Time} | 0.093 | 0.83 | 0.96 | 107 | consistent with DAG |
| Benthos -> N | Benthos_lag | N_lag | {Time} | 0.110 | 0.74 | 0.95 | 107 | consistent with DAG |
| Benthos -> N | Benthos_lag | OtherFish | {Benthos, Time} | 0.161 | 0.4 | 0.75 | 114 | consistent with DAG |
| Benthos -> N | Benthos_lag | Pmax | {Benthos, DayTemp, Flow, N} | 0.687 | 0.035 | 0.3 | 17 | consistent with DAG |
| Benthos -> N | Benthos_lag | Pmax | {Benthos, DayTemp, N, Season} | 0.687 | 0.035 | 0.3 | 17 | consistent with DAG |
| Benthos -> N | Benthos_lag | Pmax | {Benthos, DHW, N} | 0.546 | 0.19 | 0.56 | 17 | consistent with DAG |
| Benthos -> N | Benthos_lag | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.719 | 0.015 | 0.18 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Pmax | {Benthos, DayTemp, Herb, Season, Time} | 0.727 | 0.013 | 0.18 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Pmax | {Benthos, DHW, Herb, Time} | 0.574 | 0.12 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.655 | 0.044 | 0.34 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.665 | 0.038 | 0.31 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Rd | {Benthos, DHW, Herb, OtherFish} | 0.686 | 0.027 | 0.27 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.723 | 0.013 | 0.18 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.732 | 0.011 | 0.18 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Rd | {Benthos, DHW, Herb, Time} | 0.738 | 0.0099 | 0.18 | 18 | consistent with DAG |
| Benthos -> N | Benthos_lag | Season | {} | 0.357 | 0.58 | 0.89 | 18 | consistent with DAG |
| Benthos -> N | COTS | Corall | {Benthos} | 0.157 | 0.087 | 0.42 | 120 | consistent with DAG |
| Benthos -> N | COTS | Cyclone | {} | 0.002 | 0.98 | 0.99 | 120 | consistent with DAG |
| Benthos -> N | COTS | DHW | {Time} | -0.091 | 0.32 | 0.65 | 120 | consistent with DAG |
| Benthos -> N | COTS | DayTemp | {DHW} | 0.287 | 0.25 | 0.59 | 18 | consistent with DAG |
| Benthos -> N | COTS | DayTemp | {Time} | 0.322 | 0.19 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | COTS | Flow | {} | 0.319 | 0.2 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | COTS | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | 0.093 | 0.33 | 0.65 | 114 | consistent with DAG |
| Benthos -> N | COTS | Herb_lag | {Time} | -0.102 | 0.28 | 0.63 | 114 | consistent with DAG |
| Benthos -> N | COTS | N | {Benthos, Herb, Time} | 0.303 | 0.0015 | 0.046 | 107 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | COTS | N | {Benthos, Benthos_lag, Herb_lag, Time} | 0.341 | 0.00033 | 0.018 | 107 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | COTS | N_lag | {Time} | 0.495 | 5.9e-08 | 1.3e-05 | 107 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | COTS | OtherFish | {Benthos, Time} | 0.018 | 0.85 | 0.96 | 120 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, DayTemp, Flow, N} | -0.411 | 0.1 | 0.43 | 17 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, DayTemp, N, Season} | -0.369 | 0.15 | 0.49 | 17 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, DHW, N} | -0.581 | 0.014 | 0.18 | 17 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | -0.456 | 0.057 | 0.36 | 18 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, DayTemp, Herb, Season, Time} | -0.369 | 0.13 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, DHW, Herb, Time} | -0.568 | 0.014 | 0.18 | 18 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.484 | 0.042 | 0.33 | 18 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.470 | 0.049 | 0.35 | 18 | consistent with DAG |
| Benthos -> N | COTS | Pmax | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.703 | 0.0011 | 0.04 | 18 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | COTS | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.137 | 0.59 | 0.89 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.079 | 0.75 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, DHW, Herb, OtherFish} | 0.197 | 0.43 | 0.78 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.054 | 0.83 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.041 | 0.87 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, DHW, Herb, Time} | 0.188 | 0.45 | 0.79 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | 0.038 | 0.88 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | 0.050 | 0.84 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | COTS | Rd | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.168 | 0.51 | 0.85 | 18 | consistent with DAG |
| Benthos -> N | COTS | Season | {} | 0.367 | 0.13 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | Corall | Cyclone | {Benthos} | -0.093 | 0.31 | 0.65 | 120 | consistent with DAG |
| Benthos -> N | Corall | DHW | {Benthos} | 0.000 | 1 | 1 | 120 | consistent with DAG |
| Benthos -> N | Corall | DayTemp | {DHW} | 0.049 | 0.85 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Corall | DayTemp | {Benthos} | 0.348 | 0.16 | 0.5 | 18 | consistent with DAG |
| Benthos -> N | Corall | Flow | {} | 0.021 | 0.94 | 0.98 | 18 | consistent with DAG |
| Benthos -> N | Corall | Herb | {Benthos} | 0.113 | 0.22 | 0.57 | 120 | consistent with DAG |
| Benthos -> N | Corall | Herb_lag | {Benthos} | -0.112 | 0.23 | 0.59 | 114 | consistent with DAG |
| Benthos -> N | Corall | N_lag | {Benthos} | 0.022 | 0.82 | 0.96 | 107 | consistent with DAG |
| Benthos -> N | Corall | OtherFish | {Benthos} | 0.157 | 0.087 | 0.42 | 120 | consistent with DAG |
| Benthos -> N | Corall | Pmax | {Benthos, DayTemp, Flow, N} | 0.125 | 0.63 | 0.93 | 17 | consistent with DAG |
| Benthos -> N | Corall | Pmax | {Benthos, DayTemp, N, Season} | 0.127 | 0.63 | 0.93 | 17 | consistent with DAG |
| Benthos -> N | Corall | Pmax | {Benthos, DHW, N} | -0.041 | 0.88 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Corall | Pmax | {Benthos, Benthos_lag, Herb_lag, N, Time} | 0.058 | 0.82 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Corall | Pmax | {Benthos, Herb, N, Time} | 0.035 | 0.89 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Corall | Season | {} | -0.118 | 0.64 | 0.93 | 18 | consistent with DAG |
| Benthos -> N | Corall | Time | {Benthos} | -0.134 | 0.14 | 0.49 | 120 | consistent with DAG |
| Benthos -> N | Cyclone | DHW | {} | -0.110 | 0.23 | 0.59 | 120 | consistent with DAG |
| Benthos -> N | Cyclone | DayTemp | {} | -0.630 | 0.005 | 0.13 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Flow | {} | -0.003 | 0.99 | 0.99 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | -0.073 | 0.44 | 0.79 | 114 | consistent with DAG |
| Benthos -> N | Cyclone | Herb_lag | {} | -0.179 | 0.057 | 0.36 | 114 | consistent with DAG |
| Benthos -> N | Cyclone | N | {Benthos, Herb, Time} | 0.012 | 0.9 | 0.96 | 107 | consistent with DAG |
| Benthos -> N | Cyclone | N | {Benthos, Benthos_lag, Herb_lag, Time} | 0.045 | 0.64 | 0.93 | 107 | consistent with DAG |
| Benthos -> N | Cyclone | N_lag | {} | 0.051 | 0.6 | 0.91 | 107 | consistent with DAG |
| Benthos -> N | Cyclone | OtherFish | {Benthos, Time} | -0.210 | 0.021 | 0.24 | 120 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, DayTemp, Flow, N} | -0.053 | 0.84 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, DayTemp, N, Season} | -0.067 | 0.8 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.004 | 0.99 | 0.99 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, DayTemp, Herb, Season, Time} | -0.019 | 0.94 | 0.98 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.053 | 0.83 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.140 | 0.58 | 0.89 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, DHW, N} | 0.366 | 0.15 | 0.49 | 17 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, DHW, Herb, Time} | 0.384 | 0.12 | 0.48 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Pmax | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | 0.402 | 0.098 | 0.43 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | -0.354 | 0.15 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | -0.379 | 0.12 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, DayTemp, Flow, Herb, Time} | -0.271 | 0.28 | 0.63 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, DayTemp, Herb, Season, Time} | -0.267 | 0.28 | 0.63 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.406 | 0.095 | 0.42 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.429 | 0.075 | 0.42 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, DHW, Herb, OtherFish} | -0.432 | 0.074 | 0.42 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, DHW, Herb, Time} | -0.356 | 0.15 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Rd | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.554 | 0.017 | 0.2 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Season | {} | -0.367 | 0.13 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | Cyclone | Time | {} | -0.179 | 0.05 | 0.35 | 120 | consistent with DAG |
| Benthos -> N | DHW | Flow | {} | 0.045 | 0.86 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | DHW | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | -0.006 | 0.95 | 0.99 | 114 | consistent with DAG |
| Benthos -> N | DHW | Herb_lag | {Time} | -0.014 | 0.89 | 0.96 | 114 | consistent with DAG |
| Benthos -> N | DHW | N | {Benthos, Herb, Time} | 0.164 | 0.091 | 0.42 | 107 | consistent with DAG |
| Benthos -> N | DHW | N | {Benthos, Benthos_lag, Herb_lag, Time} | 0.184 | 0.057 | 0.36 | 107 | consistent with DAG |
| Benthos -> N | DHW | N_lag | {Time} | -0.026 | 0.79 | 0.96 | 107 | consistent with DAG |
| Benthos -> N | DHW | OtherFish | {Benthos, Time} | 0.031 | 0.74 | 0.95 | 120 | consistent with DAG |
| Benthos -> N | DHW | Pmax | {Benthos, DayTemp, Flow, N} | 0.134 | 0.61 | 0.91 | 17 | consistent with DAG |
| Benthos -> N | DHW | Pmax | {Benthos, DayTemp, N, Season} | 0.255 | 0.32 | 0.65 | 17 | consistent with DAG |
| Benthos -> N | DHW | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.056 | 0.83 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | DHW | Pmax | {Benthos, DayTemp, Herb, Season, Time} | 0.171 | 0.5 | 0.85 | 18 | consistent with DAG |
| Benthos -> N | DHW | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.332 | 0.18 | 0.54 | 18 | consistent with DAG |
| Benthos -> N | DHW | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.295 | 0.23 | 0.59 | 18 | consistent with DAG |
| Benthos -> N | DHW | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.242 | 0.33 | 0.65 | 18 | consistent with DAG |
| Benthos -> N | DHW | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.314 | 0.2 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | DHW | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.269 | 0.28 | 0.63 | 18 | consistent with DAG |
| Benthos -> N | DHW | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.321 | 0.19 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | DHW | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.285 | 0.25 | 0.59 | 18 | consistent with DAG |
| Benthos -> N | DHW | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.306 | 0.22 | 0.57 | 18 | consistent with DAG |
| Benthos -> N | DHW | Season | {} | 0.502 | 0.034 | 0.3 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | Flow | {Season} | 0.097 | 0.7 | 0.95 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | Herb | {Benthos, Benthos_lag, Herb_lag, Time} | 0.202 | 0.42 | 0.77 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | Herb | {DHW} | 0.282 | 0.26 | 0.6 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | Herb_lag | {Time} | 0.423 | 0.08 | 0.42 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | Herb_lag | {DHW} | 0.359 | 0.14 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | N | {Benthos, Herb, Time} | 0.046 | 0.86 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | DayTemp | N | {Benthos, Benthos_lag, Herb_lag, Time} | 0.166 | 0.52 | 0.86 | 17 | consistent with DAG |
| Benthos -> N | DayTemp | N | {DHW} | -0.087 | 0.74 | 0.95 | 17 | consistent with DAG |
| Benthos -> N | DayTemp | N_lag | {Time} | 0.358 | 0.14 | 0.49 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | N_lag | {DHW} | 0.030 | 0.9 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | OtherFish | {Benthos, Time} | 0.240 | 0.34 | 0.65 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | OtherFish | {DHW} | 0.254 | 0.31 | 0.65 | 18 | consistent with DAG |
| Benthos -> N | DayTemp | Time | {DHW} | 0.137 | 0.59 | 0.89 | 18 | consistent with DAG |
| Benthos -> N | Flow | Herb | {} | 0.213 | 0.4 | 0.74 | 18 | consistent with DAG |
| Benthos -> N | Flow | Herb_lag | {} | 0.274 | 0.27 | 0.62 | 18 | consistent with DAG |
| Benthos -> N | Flow | N | {} | -0.118 | 0.65 | 0.93 | 17 | consistent with DAG |
| Benthos -> N | Flow | N_lag | {} | -0.073 | 0.77 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Flow | OtherFish | {} | -0.230 | 0.36 | 0.68 | 18 | consistent with DAG |
| Benthos -> N | Flow | Time | {} | 0.105 | 0.68 | 0.94 | 18 | consistent with DAG |
| Benthos -> N | Herb | N_lag | {Benthos, Benthos_lag, Herb_lag, Time} | -0.114 | 0.24 | 0.59 | 107 | consistent with DAG |
| Benthos -> N | Herb | OtherFish | {Benthos, Time} | 0.249 | 0.006 | 0.13 | 120 | consistent with DAG |
| Benthos -> N | Herb | Pmax | {Benthos, DayTemp, Flow, N} | 0.066 | 0.8 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Herb | Pmax | {Benthos, DayTemp, N, Season} | 0.112 | 0.67 | 0.93 | 17 | consistent with DAG |
| Benthos -> N | Herb | Pmax | {Benthos, DHW, N} | -0.154 | 0.56 | 0.89 | 17 | consistent with DAG |
| Benthos -> N | Herb | Pmax | {Benthos, Benthos_lag, Herb_lag, N, Time} | -0.005 | 0.98 | 0.99 | 17 | consistent with DAG |
| Benthos -> N | Herb | Season | {} | 0.410 | 0.091 | 0.42 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | N | {Benthos, Herb, Time} | 0.267 | 0.0054 | 0.13 | 107 | consistent with DAG |
| Benthos -> N | Herb_lag | N_lag | {Time} | 0.017 | 0.86 | 0.96 | 107 | consistent with DAG |
| Benthos -> N | Herb_lag | OtherFish | {Benthos, Time} | -0.011 | 0.9 | 0.96 | 114 | consistent with DAG |
| Benthos -> N | Herb_lag | Pmax | {Benthos, DayTemp, Flow, N} | 0.382 | 0.13 | 0.49 | 17 | consistent with DAG |
| Benthos -> N | Herb_lag | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.292 | 0.24 | 0.59 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Pmax | {Benthos, DayTemp, N, Season} | 0.368 | 0.15 | 0.49 | 17 | consistent with DAG |
| Benthos -> N | Herb_lag | Pmax | {Benthos, DayTemp, Herb, Season, Time} | 0.256 | 0.3 | 0.65 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Pmax | {Benthos, DHW, N} | 0.194 | 0.45 | 0.79 | 17 | consistent with DAG |
| Benthos -> N | Herb_lag | Pmax | {Benthos, DHW, Herb, Time} | 0.159 | 0.53 | 0.86 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.228 | 0.36 | 0.68 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Rd | {Benthos, DayTemp, Flow, Herb, Time} | 0.287 | 0.25 | 0.59 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | 0.240 | 0.34 | 0.65 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Rd | {Benthos, DayTemp, Herb, Season, Time} | 0.312 | 0.21 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Rd | {Benthos, DHW, Herb, OtherFish} | 0.347 | 0.16 | 0.5 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Rd | {Benthos, DHW, Herb, Time} | 0.414 | 0.087 | 0.42 | 18 | consistent with DAG |
| Benthos -> N | Herb_lag | Season | {} | 0.233 | 0.35 | 0.67 | 18 | consistent with DAG |
| Benthos -> N | N | N_lag | {Benthos, Benthos_lag, Herb_lag, Time} | 0.037 | 0.71 | 0.95 | 101 | consistent with DAG |
| Benthos -> N | N | N_lag | {Benthos, Herb, Time} | 0.015 | 0.88 | 0.96 | 101 | consistent with DAG |
| Benthos -> N | N | Rd | {Benthos, Corall, DayTemp, Flow, Herb, OtherFish} | 0.248 | 0.34 | 0.65 | 17 | consistent with DAG |
| Benthos -> N | N | Rd | {Benthos, Corall, DayTemp, Herb, OtherFish, Season} | 0.255 | 0.32 | 0.65 | 17 | consistent with DAG |
| Benthos -> N | N | Rd | {Benthos, Corall, DHW, Herb, OtherFish} | 0.193 | 0.46 | 0.79 | 17 | consistent with DAG |
| Benthos -> N | N | Rd | {Benthos, Corall, Herb, OtherFish, Time} | 0.526 | 0.03 | 0.29 | 17 | consistent with DAG |
| Benthos -> N | N | Season | {} | -0.158 | 0.55 | 0.88 | 17 | consistent with DAG |
| Benthos -> N | N_lag | OtherFish | {Benthos, Time} | 0.220 | 0.023 | 0.24 | 107 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, DayTemp, Flow, N} | -0.015 | 0.95 | 0.99 | 17 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, DayTemp, Flow, Herb, Time} | 0.004 | 0.99 | 0.99 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | 0.093 | 0.71 | 0.95 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, DayTemp, N, Season} | -0.123 | 0.64 | 0.93 | 17 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, DayTemp, Herb, Season, Time} | -0.141 | 0.58 | 0.89 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | 0.007 | 0.98 | 0.99 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, DHW, N} | -0.117 | 0.65 | 0.93 | 17 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, DHW, Herb, Time} | -0.159 | 0.53 | 0.86 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Pmax | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.166 | 0.51 | 0.85 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, DayTemp, Flow, Herb, OtherFish} | -0.070 | 0.78 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, DayTemp, Flow, Herb, Time} | -0.146 | 0.56 | 0.89 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.256 | 0.3 | 0.65 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, DayTemp, Herb, OtherFish, Season} | -0.101 | 0.69 | 0.95 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, DayTemp, Herb, Season, Time} | -0.175 | 0.49 | 0.83 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, Benthos_lag, DayTemp, Herb_lag, Season, Time} | -0.341 | 0.17 | 0.52 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, DHW, Herb, OtherFish} | 0.094 | 0.71 | 0.95 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, DHW, Herb, Time} | 0.058 | 0.82 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Rd | {Benthos, Benthos_lag, DHW, Herb_lag, Time} | -0.341 | 0.17 | 0.52 | 18 | consistent with DAG |
| Benthos -> N | N_lag | Season | {} | -0.324 | 0.19 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | OtherFish | Pmax | {Benthos, DayTemp, Flow, N} | -0.444 | 0.074 | 0.42 | 17 | consistent with DAG |
| Benthos -> N | OtherFish | Pmax | {Benthos, DayTemp, N, Season} | -0.432 | 0.083 | 0.42 | 17 | consistent with DAG |
| Benthos -> N | OtherFish | Pmax | {Benthos, DHW, N} | -0.422 | 0.091 | 0.42 | 17 | consistent with DAG |
| Benthos -> N | OtherFish | Pmax | {Benthos, Benthos_lag, Herb_lag, N, Time} | -0.453 | 0.068 | 0.41 | 17 | consistent with DAG |
| Benthos -> N | OtherFish | Pmax | {Benthos, Herb, N, Time} | -0.490 | 0.046 | 0.34 | 17 | consistent with DAG |
| Benthos -> N | OtherFish | Season | {} | 0.096 | 0.7 | 0.95 | 18 | consistent with DAG |
| Benthos -> N | Pmax | Rd | {Benthos, Corall, DayTemp, Flow, Herb, OtherFish} | 0.762 | 0.00024 | 0.017 | 18 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | Pmax | Rd | {Benthos, DayTemp, Flow, N} | 0.777 | 0.00024 | 0.017 | 17 | VIOLATED (BH-adjusted p < 0.05) |
| Benthos -> N | Pmax | Season | {DHW, DayTemp, Flow} | -0.313 | 0.21 | 0.56 | 18 | consistent with DAG |
| Benthos -> N | Pmax | Season | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.111 | 0.66 | 0.93 | 18 | consistent with DAG |
| Benthos -> N | Pmax | Season | {Benthos, DayTemp, Flow, Herb, Time} | -0.187 | 0.46 | 0.79 | 18 | consistent with DAG |
| Benthos -> N | Pmax | Season | {Benthos, DayTemp, Flow, N} | -0.146 | 0.58 | 0.89 | 17 | consistent with DAG |
| Benthos -> N | Pmax | Time | {Benthos, DHW, N} | -0.092 | 0.73 | 0.95 | 17 | consistent with DAG |
| Benthos -> N | Pmax | Time | {Benthos, DayTemp, N, Season} | 0.096 | 0.71 | 0.95 | 17 | consistent with DAG |
| Benthos -> N | Pmax | Time | {Benthos, DayTemp, Flow, N} | 0.074 | 0.78 | 0.96 | 17 | consistent with DAG |
| Benthos -> N | Rd | Season | {DHW, DayTemp, Flow} | -0.081 | 0.75 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Rd | Season | {Benthos, Benthos_lag, DayTemp, Flow, Herb_lag, Time} | -0.088 | 0.73 | 0.95 | 18 | consistent with DAG |
| Benthos -> N | Rd | Season | {Benthos, DayTemp, Flow, Herb, Time} | 0.021 | 0.93 | 0.98 | 18 | consistent with DAG |
| Benthos -> N | Rd | Season | {Benthos, DayTemp, Flow, Herb, OtherFish} | 0.034 | 0.89 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Rd | Time | {Benthos, DHW, Herb, OtherFish} | -0.042 | 0.87 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Rd | Time | {Benthos, DayTemp, Herb, OtherFish, Season} | -0.025 | 0.92 | 0.98 | 18 | consistent with DAG |
| Benthos -> N | Rd | Time | {Benthos, DayTemp, Flow, Herb, OtherFish} | -0.051 | 0.84 | 0.96 | 18 | consistent with DAG |
| Benthos -> N | Season | Time | {} | 0.314 | 0.2 | 0.56 | 18 | consistent with DAG |
