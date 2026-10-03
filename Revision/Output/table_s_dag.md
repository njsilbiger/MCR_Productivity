# Table S-DAG: minimal adjustment sets for every causal estimand

Two DAG variants are carried through this table pending the Task 4.3
d-separation test against data (not yet run -- see note at the end of
`03_dag.R`): `Rd -> N` (respiration-linked remineralisation) vs
`Benthos -> N` (nitrogen driven directly by standing benthic biomass).
Decided with the PI 2026-10-01 to draw and tabulate both rather than
choose one now.

Each row is one causal estimand (Task 4.2); each gets its own model in
Section 7 using only its own adjustment set (Arif & MacNeil 2023) -- an
estimand's effect is never read off a model fit for a different estimand.
E11 (`Benthos -> Corallivores / Herbivores`) is split into E11a/E11b
because dagitty's adjustment sets are defined per exposure-outcome pair.

| DAG variant | ID | Estimand | Effect | Adjustment set(s) | # minimal sets |
|---|---|---|---|---|---|
| Rd -> N | E1 | DHW -> Benthos (coral share) | total | {Time} | 1 |
| Rd -> N | E2 | COTS -> Benthos | total | {Benthos_lag, Time} | 1 |
| Rd -> N | E3 | Herb_lag -> Benthos (recovery mediation) | total | {Time} | 1 |
| Rd -> N | E4 | Benthos -> Rd (direct, not via fish) | direct | {Corall, DayTemp, Flow, Herb, OtherFish}; {Corall, DayTemp, Herb, OtherFish, Season}; {Corall, DHW, Herb, OtherFish} | 3 |
| Rd -> N | E5 | Benthos -> Rd (total) | total | {Benthos_lag, DayTemp, Flow, Herb_lag, Time}; {Benthos_lag, DayTemp, Herb_lag, Season, Time}; {Benthos_lag, DHW, Herb_lag, Time} | 3 |
| Rd -> N | E6 | Fish (Herb+Corall+OtherFish) -> Rd (direct) | direct | {Benthos, DayTemp, Flow}; {Benthos, DayTemp, Season}; {Benthos, DHW}; {Benthos, Benthos_lag, Herb_lag, Time} | 4 |
| Rd -> N | E7 | DayTemp -> Rd (physiological, direct) | direct | {Benthos, Flow, Herb, OtherFish}; {Benthos, Herb, OtherFish, Season}; {Benthos, Flow, Herb, Time}; {Benthos, Herb, Season, Time}; {Benthos, Benthos_lag, Flow, Herb_lag, Time}; {Benthos, Benthos_lag, Herb_lag, Season, Time}; {DHW, Flow}; {DHW, Season} | 8 |
| Rd -> N | E8 | DHW -> Rd (total) | total | {Time} | 1 |
| Rd -> N | E9 | Benthos -> Pmax (direct) | direct | {DayTemp, Flow, N}; {Corall, DayTemp, Flow, Herb, OtherFish, Rd, Time} | 2 |
| Rd -> N | E10 | DayTemp -> Pmax (direct) | direct | {Benthos, Flow, N}; {Benthos, Corall, Flow, Herb, OtherFish, Rd, Time}; {Benthos, Corall, DHW, Flow, Herb, OtherFish, Rd} | 3 |
| Rd -> N | E11a | Benthos -> Herbivores (total) | total | {Benthos_lag, Herb_lag, Time} | 1 |
| Rd -> N | E11b | Benthos -> Corallivores (total) | total | {} | 1 |
| Rd -> N | E12 | Fish (Herb+Corall+OtherFish) -> N (total) | total | {Benthos, DayTemp, Flow, Time}; {Benthos, DayTemp, Season, Time}; {Benthos, DHW, Time}; {Benthos, Benthos_lag, Herb_lag, Time} | 4 |
| Benthos -> N | E1 | DHW -> Benthos (coral share) | total | {Time} | 1 |
| Benthos -> N | E2 | COTS -> Benthos | total | {Benthos_lag, Time} | 1 |
| Benthos -> N | E3 | Herb_lag -> Benthos (recovery mediation) | total | {Time} | 1 |
| Benthos -> N | E4 | Benthos -> Rd (direct, not via fish) | direct | {Corall, DayTemp, Flow, Herb, OtherFish}; {Corall, DayTemp, Herb, OtherFish, Season}; {Corall, DHW, Herb, OtherFish} | 3 |
| Benthos -> N | E5 | Benthos -> Rd (total) | total | {Benthos_lag, DayTemp, Flow, Herb_lag, Time}; {Benthos_lag, DayTemp, Herb_lag, Season, Time}; {Benthos_lag, DHW, Herb_lag, Time} | 3 |
| Benthos -> N | E6 | Fish (Herb+Corall+OtherFish) -> Rd (direct) | direct | {Benthos, DayTemp, Flow}; {Benthos, DayTemp, Season}; {Benthos, DHW}; {Benthos, Benthos_lag, Herb_lag, Time} | 4 |
| Benthos -> N | E7 | DayTemp -> Rd (physiological, direct) | direct | {Benthos, Flow, Herb, OtherFish}; {Benthos, Herb, OtherFish, Season}; {Benthos, Flow, Herb, Time}; {Benthos, Herb, Season, Time}; {Benthos, Benthos_lag, Flow, Herb_lag, Time}; {Benthos, Benthos_lag, Herb_lag, Season, Time}; {DHW, Flow}; {DHW, Season} | 8 |
| Benthos -> N | E8 | DHW -> Rd (total) | total | {Time} | 1 |
| Benthos -> N | E9 | Benthos -> Pmax (direct) | direct | {DayTemp, Flow, N}; {DayTemp, N, Season}; {DHW, N} | 3 |
| Benthos -> N | E10 | DayTemp -> Pmax (direct) | direct | {Benthos, Flow, N}; {Benthos, N, Season}; {Benthos, Flow, Herb, Time}; {Benthos, Herb, Season, Time}; {Benthos, Benthos_lag, Flow, Herb_lag, Time}; {Benthos, Benthos_lag, Herb_lag, Season, Time}; {DHW, Flow}; {DHW, Season} | 8 |
| Benthos -> N | E11a | Benthos -> Herbivores (total) | total | {Benthos_lag, Herb_lag, Time} | 1 |
| Benthos -> N | E11b | Benthos -> Corallivores (total) | total | {} | 1 |
| Benthos -> N | E12 | Fish (Herb+Corall+OtherFish) -> N (total) | total | {Benthos, Time} | 1 |
