[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Hydrogen Network Expansion Planning considering the Chicken-and-egg Dilemma and Market Uncertainty

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Hydrogen Network Expansion Planning considering the Chicken-and-egg Dilemma and Market Uncertainty](https://doi.org/10.1287/ijoc.2025.1111) by Beste Basciftci, Sezen Ece Kayacık, Gerben B. Ouwersloot, Albert H. Schrotenboer, Iris F. A. Vis, and Evrim Ursavas.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2025.1111

https://doi.org/10.1287/ijoc.2025.1111.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{CacheTest,
  author =        {Basciftci, Beste and Kayacık, Sezen Ece and Ouwersloot, Gerben B. and Schrotenboer, Albert H. and Vis, Iris F. A. and Ursavas, Evrim},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Hydrogen Network Expansion Planning considering the Chicken-and-egg Dilemma and Market Uncertainty}},
  year =          {2025},
  doi =           {10.1287/ijoc.2025.1111.cd},
  url =           {https://github.com/INFORMSJoC/2025.1111},
  note =          {Available for download at https://github.com/INFORMSJoC/2025.1111},
}  
```

## Description

This software runs all numerical experiments, comparing different setups and enhancements of the column-and-constraint generation algorithm. It supports the figures and tables found in section 4.

## Building

The code is written in Python (3.12.6), using the Gurobi (12.0.3) solver. Its dependencies are numpy (2.3.2), gurobipy (13.0.1), scipy (1.16.1) and joblib (1.5.2).


## Results

Figure 1 and Tables 1, 2, and 3 show performance metrics under different setups, and are based on the data in results/data.csv.

Table 4 shows performance metrics under different budget uncertainty sets, and are based on the data in results/budget_set_data.csv.

Figure 3 shows the production capacity planning of the case study with and without the consideration of DDU.

Figure 4 shows the effect of different DDU levels on the production capacity planning.

Figure 5 shows the effect of different budget uncertainty sets on the production capacity planning.

Tables 5 and 6 show-of-sample results of the case study under the normal and beta distributions.

Figure 6 shows how deviations from the expected DDU level affects profits under different scenario simulations.

Figure 7 shows the production capacity planning of the case study under a lower import cost with and without the consideration of DDU.

The data for Figures 3, 4, 5, and 7 can be found in results/case_study_investments.csv.

The data for Tables 5 and 6, and Figure 6 can be found in results/case_study_simulation.csv.

## Replicating

To reproduce a row of the main computation results in results/data.csv, do

```
python computational.py --job {JOB_NUMBER}
```

Where {JOB_NUMBER} is the index of the job in the list created in the script as JOBlist.

Similarly, to reproduce the data of budget_set_data.csv, run

```
python computational_budget.py --job {JOB_NUMBER}
```

For all runs in  case_study_investments.csv or case_study_simulation_data.csv with a location-based first moment function, run

```
python simulation.py --job {JOB_NUMBER}
```

And lastly, for all runs in case_study_investments.csv with a capacity-based first moment function, run

```
python simulation_cap_ddu.py --job {JOB_NUMBER}
```
