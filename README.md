# Powering-tropical-cities

This repository contains the code and minimal data necessary to reproduce the key computational analyses presented in the manuscript **“Powering tropical cities: Photovoltaics at scale via decentralized e-mobility charging”** (submitted to *Nature Communications*).

All materials are provided in anonymized form for peer-review purposes only.

---

## 1. Repository Structure

```text
Powering-Tropical-Cities/
├── code/                          # All analysis and simulation scripts
│   ├── Figure1_4Sharing.m           # Figure 1 plotting, the data used is included in the code; 
│   ├── Figure2_4Sharing.m           # Figure 2 plotting, the data used is included in the code; Also, Fig S13, S14 are included
│   ├── Figure3_4Sharing.m           # Figure 3 plotting, the data used is included in the code
│   ├── Figure4_4Sharing.m           # Figure 4 plotting, the data used is included in the code
│   ├── Figure4_4Sharing.m           # Figure 4 plotting, the data used is included in the code
│   ├── Figure4_4Sharing.m           # Figure 4 plotting, the data used is included in the code
│   ├── DC_PowerFlow_4Sharing.m          # DC power flow simulation
│   ├── AC_PowerFlow_4Sharing.m          # AC power flow simulation
│   ├── ElectricityGridModel_n_Map_4Sharing.m          # Electricity Grid Model visualisation (DC power flow modeling)
│   ├── xxxx          # xxxx
│   ├── xxxx          # xxx
│   └── Helper.m          # xxxx
│       ├── load_data.m
│       ├── save_results.m
│       └── plotting_utils.m
│
├── data/                          # Input datasets (synthetic/anonymized)
│   ├── Data4Figure           # Partial data used to repproduce the figures
│   ├── LineLoad              # Line loading data
│   ├── ParkedVehNum          # Parked Veh Num data (for Main Fig. 4)
│   ├── PowerGrid             # Electricity Grid Topology
│   ├── SystemLoad            # System demand (city-wide demand)
│   ├── mp14-plng-area-no-sea-planning-area         # planning area shp file, available at: https://data.gov.sg/collections/1700/view
│   └── xxxxxxx               # xxxxx
│
├── results/                       # Outputs generated after running all "*.m"
│   ├── figure_1...
│   ├── figure_2...
│   ├── figure_3...
│   ├── figure_4...
│   └── ...
│
└── README.md                      # This file
```


---

## 2. How to Use

### Step 1 — Set up environment
The analyses were implemented primarily in **MATLAB R2025a**.  

### Step 2 — Run workflow
run('*.m')

All generated figures and summary tables will appear in the results/ folder.

## Notes:
i). The code for **Battery degradation simulation** is adapted from BLAST-List (https://github.com/NREL/BLAST-Lite). All the revelant battery models and parameters can be found in that repository. 

ii). 

