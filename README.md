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
│   ├── EV_district_level_optimization_4Sharing.m          # EV district-level optimization demonstration model
│   ├── EV_system_level_optimization_4Sharing.m          # EV system-level optimization demonstration model
│
├── data/                          # Input datasets (synthetic/anonymized)
│   ├── Data4Figure           # Main data used to reproduce the figures
│   ├── LineLoad              # Line loading data
│   ├── ParkedVehNum          # Parked Veh Num data (for Main Fig. 4)
│   ├── PowerGrid             # Electricity network Topology (and parameters)
│   ├── SystemLoad            # System demand (city-wide demand)
│   ├── mp14-plng-area-no-sea-planning-area         # planning area shp file, available at: https://data.gov.sg/collections/1700/view
│   └── MobilityInputSample_dummy               # Dummy input for the EV charging optimization model
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



## 2. Requirements

This AC power flow modeling code is implemented in MATLAB and relies on the MATPOWER toolbox.

### MATLAB
- MATLAB (R2023b or later recommended)  

### MATPOWER
- MATPOWER power system simulation package (a free, open-source tools for electric power system simulation and optimization)
- Tested with *MATPOWER 8.1* 
  (<https://matpower.org>)
- MATPOWER must be installed and added to the MATLAB path.
  For example:

```matlab
addpath('path_to_matpower');
savepath;
```


## 2. How to Use

### Step 1 — Set up environment


### Step 2 — Run workflow
run('*.m')

All generated figures and summary tables will appear in the *results* folder.

## Notes:
i). The code for **Battery degradation simulation** is adapted from BLAST-List (https://github.com/NREL/BLAST-Lite). All the revelant battery models and parameters can be found in that repository. 

ii). Electricity demand and solar irradiance data are publicly available from the Energy Market Authority of Singapore (https://www.ema.gov.sg/resources/statistics). Data for the derivation of location-specific PV generation potentials are also publicly available as described in [1]. Travel survey and census data are publicly available from the Singapore Department of Statistics (https://www.singstat.gov.sg/publications/reference/cop2020). Individual-level mobility data are not publicly available due to privacy considerations, while the data to reproduce the findings and figures in this paper are all available in this repository. 

[1]. Caviezel, D., Waibel, C., Schläpfer, M. & Schlueter, A. Vehicle-to-grid coupled photo-voltaic optimization for Singapore at a district resolution. In 36th International Conference on Efficiency, Cost, Optimization, Simulation and Environmental Impact of Energy Systems (ECOS), 3327–3338 (2023).574

