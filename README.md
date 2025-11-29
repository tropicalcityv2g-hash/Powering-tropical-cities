# Powering-tropical-cities

This repository contains the code and minimal data necessary to reproduce the key computational analyses presented in the manuscript **“Powering tropical cities: Photovoltaics at scale via decentralized e-mobility charging”** (submitted to *Nature Communications*).

All materials are provided in anonymized form for peer-review purposes only.

---

## 1. Repository Structure

```text
Powering-Tropical-Cities/
├── code/                          # All analysis and simulation scripts
│   ├── main_analysis.m            # Main entry point that runs the complete workflow
│   ├── optimization_model.m       # Decentralized EV–PV–grid optimization model
│   ├── powerflow_simulation.m     # DC power-flow computation for each scenario
│   ├── plot_generation.m          # Figure creation scripts for main text and SI
│   └── helper_functions/          # Utility functions used by other scripts
│       ├── load_data.m
│       ├── save_results.m
│       └── plotting_utils.m
│
├── data/                          # Input datasets (synthetic/anonymized)
│   ├── load_profiles_sample.mat           # Baseline demand profiles (hourly)
│   ├── pv_scenarios_sample.mat            # PV generation profiles (BAS/ACC/MAX)
│   ├── network_parameters_sample.mat      # Simplified 66 kV grid structure
│   └── ev_adoption_sample.mat             # EV adoption and capacity data
│
├── results/                       # Outputs generated after running main_analysis.m
│   ├── figure_1_PV_matching.png
│   ├── figure_2_EV_distribution.png
│   ├── figure_3_powerflow.png
│   ├── figure_4_ramp_reduction.png
│   ├── summary_metrics.csv
│   └── intermediate_results.mat
│
├── requirements.txt               # Python dependencies (if applicable)
└── README.md                      # This file
```


---

## 2. How to Use

### Step 1 — Set up environment
The analyses were implemented primarily in **MATLAB R2025a**.  

### Step 2 — Run workflow
run('code/main_analysis.m')

All generated figures and summary tables will appear in the results/ folder.
