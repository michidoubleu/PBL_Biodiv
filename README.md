# PBL_Biodiv Post-Processing Pipeline

This repository contains the post-processing workflow for the **Biodiversity link to netcdf**, as used internally at **PBL_Biodiv**.

The pipeline performs:

1. Loading scenario-specific biodiversity link output files (`*.RData`)
2. Land-cover mapping  
   - G4M → DS land-cover classes  
   - Protected-area harmonization
3. Aggregation of LC and LULUC outputs
4. Writing LC and LULUC CSV files
5. Conversion to **netCDF**

---

## Input Data and Folder Structure

### A. Network drive inputs (main data source)

The script reads all scenario data and raw biodiversity link output files from:

```
P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results/Biodiversity_Link/
```

Required files in this directory include:

- `scenario_mapping_*.csv`
- Biodiversity link outputs:  
  `biodiversity_<SCEN>_<DATE>.RData`

These `.RData` files contain lists where element `[[5]]` holds the link results used by the script.

---




