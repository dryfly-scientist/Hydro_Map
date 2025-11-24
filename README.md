# Hydro_Map

# Terrestrial Nitrogen Retention Hydrologic Monitoring Map

This repository provides two Google Earth Engine (GEE) workflows during NASA DEVELOP Fall 2025 for evaluating **Nitrogen Retention Indices (NRI)** across the lakeshore Wisconsin counties study region, and an interactive **Terrestrial Nitrogen Retention Index (TNRI)** output for project partners that incorporates water quality data, hydrographic parameters, and hydrologic connectivity.

Reference scripts include:

1. **Hydro_Map_Base.js** — Statewide NRI (no water quality inputs; baseline nitrogen retention index)
2. **Hydro_Map_PalmerFox.js** — TNRI demonstration for the Palmer–Fox HUC12 using USGS nitrate + nitrite

Project partners (USDA NRCS) can apply this workflow to any HUC12 watershed polygon with the following inputs:
- water quality data .csv,  
- a matching streamflow discharge estimate (USGS)  
- HUC12 watershed boundary of interest

---

## Application

The goal of this workflow is to map areas of likely terrestrial nitrogen retention by combining vegetation density (NDVI), terrain (DEM + slope), in-situ water quality parameters, river channel proximity to streams, and MERIT Hydro parameters; including flow accumulation, flow direction, and hydrologically adjusted elevation + river channel width.  

This allows partners to:
- identify key riparian corridors for improved soil health, runoff mitigation, stable water quality, and carbon sequstration productivity  
- evaluate nutrient retention potential across cover crop adoption zones  
- compare watersheds' NRI & TNRI   
- Expoand workflow outside of current ROI  


##  Data Sources 

Partners can locate these in the Google Earth Engine Data Catalog.

### **Raster Datasets**
Format: / Dataset / GEE ID / Purpose /

- / Cropland Data Layer (CDL) / `USDA/NASS/CDL` / Winter cover-crop detection /
- / SRTM DEM 30m / `USGS/SRTMGL1_003` / Slope, elevation, topography /
- / MERIT Hydro / `MERIT/Hydro/v1_0_1` / Hydrography; flow accumulation weighting (primary function) /

### **Vector Datasets**
Format: / Dataset / Source in GEE / Purpose /

- / USGS HUC12 Watersheds / Search “WBD” or “HUC12” in Data Catalog / Watershed boundaries /
- / NHD Flowlines / Search “NHDFlowline” / Riparian buffer + proximity weighting /
- / Partner Water-Quality Data / Upload CSV as **Table** / Nutrient concentration (mg/L), timestamps /

---

## How Do I (you) Use This Repository?

### 1. Open Google Earth Engine Code Editor  
https://code.earthengine.google.com/

### 2. Copy each `.js` file into its own script tab within GEE
- Paste **Hydro_Map_Base.js** into Script A  
- Paste **Hydro_Map_PalmerFox.js** into Script B  

### 3. Modify only what your project requires  
- Your water-quality CSV  
- Your discharge value (m³/s) from USGS NWIS  
- Any custom date ranges
- Your watershed HUC12 (Select polygon ==>     

### 4. Run  
Outputs export as GeoTIFF to Google Drive upon selection of HUC12 watershed boundary.

---

# Overview of NRI methodology

The NRI model combines:

- **NDVI_scaled**  
- **Proximity_to_Streams** (exponential decay)  
- **Lowland Index** (DEM + slope)  
- **Flow Proxy** (hydrologic downslope accumulation)  
- **30/120/240 m riparian mask** (vegetated corridors)

### Formula:

NRI = NDVI_scaled × Proximity_to_Streams × Lowland_Index × Flow_Proxy


The base NRI output is a normalized 30 m study area-wide raster that can be used as the baseline for TNRI wherein partners integrate their own in-situ water quality parameters.
These may include nitrate + nitrite, dissolved oxygen, turbidity, phosphorus, etc.

-------------------------------------------

# TNRI Extension (water quality + discharge)

This is included in Hydro_Map_PalmerFox.js. User will input their parameters into base NRI, and select any watershed. Output will reflect the framework of PalmerFox example findings. 

### Step 1 — Upload water-quality table
Format (CSV):

system:time_start,value
- 2023-03-01T00:00:00Z,1.84 (example value of nitrate + nitrite in mg/L)
- 2023-03-14T00:00:00Z,1.92
- 2023-05-01T00:00:00Z,2.10


Upload it to GEE as:  
Assets ==> NEW ==> Table upload.  
Replace placeholder asset path in the script.

### Step 2 — Provide mean discharge
Use USGS NWIS: https://waterdata.usgs.gov/
Example from Palmer-Fox: https://waterdata.usgs.gov/monitoring-location/USGS-05545750/#dataTypeId=continuous-00060-0&period=P7D&showFieldMeasurements=true

Convert from cfs to m³/s using:
m^3/s = cfs × 0.0283168


Replace the placeholder discharge value in the script.

### Step 3 — TNRI computation
The script multiplies:
- the clipped NRI for the selected watershed,  
- MERIT Hydro flow accumulation,  
- and nutrient flux computed from concentration × discharge.

Formula:
TNRI = (NRI) * (MERIT_flow) * (nutrient_flux)


### Step 4 — Export
The script exports a TNRI GeoTIFF to Google Drive.

---

## Script Descriptions

### Hydro_Map_Base.js
- Computes Landsat NDVI (scaled SR).
- Computes slope and lowland index from SRTM.
- Creates distance-to-stream weighting using NHD Flowlines.
- Builds 30/120/240 m riparian masks.
- Produces statewide NRI and exports a GeoTIFF.

### Hydro_Map_PalmerFox.js
- Includes everything in the base script.
- Clipped example for Palmer–Fox watershed.
- Adds MERIT Hydro flow accumulation.
- Adds partner-ready nitrate + discharge module.
- Computes spatialized nitrate load and TNRI.
- Produces all required maps and exports to Drive.

---

## Reproducibility Notes

- All datasets listed are public in GEE.
- The only required partner data is a water-quality CSV and streamflow discharge from nearest USGS stream gauge.
- Discharge must match the same time window as the CSV.
- This method works for nitrogen, phosphorus, turbidity, DO, or any parameter with concentration data.

---
