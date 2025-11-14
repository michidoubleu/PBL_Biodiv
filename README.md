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

### A. Input sources

The script reads all scenario data and raw biodiversity link output files from the specified network drive:

```
P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results/Biodiversity_Link/
```

This can be adjusted to the specific folder in the main script, line 30 under `define.path`

Required files in this directory include:
- `scenario_mapping_*.csv`
- Biodiversity link outputs.
These `.RData` files contain lists where element `[[5]]` holds the link results used by the script.


Additional supporting files are provided in the template folder on the IIASA PDrive at:
```
P:/globiom/Projects/PBL_BIODIV_2025/Postprocessing_ncdf/template
```
This is not to be changed!

---

## Additional Notes

- Execution time:
  - ~10 minutes per scenario for CSV output
  - additional ~15 minutes per scenario for netCDF output  




