// **********
//    Wisconsin NRI Reference
//    Statewide NRI using NDVI (2023 reference year), DEM, slope, proximity, and 30m/120m/240m buffers.
//    Partner module (line 183+) shows how to add water quality data and discharge
//    to generate a full terrestrial nitrogen retention index (TNRI)
//    for any HUC12 watershed within study area.
// **********


// ********** Imports
var counties = table3;   // WI counties
var flow_raw = table2;   // WDNR flowlines
var huc12    = table;    // WI USGS HUC12 watersheds


// ********** ROI 
var county_list = [
  "Door","Kewaunee","Manitowoc","Sheboygan",
  "Ozaukee","Milwaukee","Racine","Kenosha"
];

var wisc_roi = counties.filter(
  ee.Filter.inList("COUNTY_NAM", county_list)
);

// Watersheds within ROI
var allWatersheds = huc12.filterBounds(wisc_roi);

// Flowlines within ROI
var flowlines = flow_raw.filterBounds(wisc_roi);


// **********
//    NDVI Calculation
// **********

function maskLandsatL2(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));
  return img.updateMask(mask);
}

function scaleSR(img) {
  var srBands = img.bandNames()
    .filter(ee.Filter.stringContains('item', 'SR_B'));
  var scaled = img.select(srBands)
    .multiply(0.0000275)
    .add(-0.2);
  return img.addBands(scaled, null, true);
}

function addNDVI(img) {
  var bands = img.bandNames();
  var hasB5 = bands.contains('SR_B5');
  var nir = ee.String(ee.Algorithms.If(hasB5, 'SR_B5', 'SR_B4'));
  var red = ee.String(ee.Algorithms.If(hasB5, 'SR_B4', 'SR_B3'));
  return img.addBands(
    img.normalizedDifference([nir, red]).rename('NDVI')
  );
}

//  NDVI (adjust date window based on partner study period)
function getSpringNDVI_2023() {
  var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
  var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");

  return l8.merge(l9)
    .filterBounds(wisc_roi)
    .filterDate('2023-02-01','2023-05-31')
    .map(maskLandsatL2)
    .map(scaleSR)
    .map(addNDVI)
    .select('NDVI')
    .median()
    .clip(wisc_roi);
}

var springNDVI_2023 = getSpringNDVI_2023();
var ndvi_scaled = springNDVI_2023.unitScale(0, 0.8).clamp(0, 1);


// **********
//    DEM + Terrain
// **********

var dem = ee.Image("USGS/SRTMGL1_003").clip(wisc_roi);
var slope = ee.Terrain.slope(dem);

var elev_scaled = dem.unitScale(100, 400).clamp(0, 1);
var slope_norm  = slope.unitScale(0, 15);

// Flow accumulation proxy from slope
var flow_proxy = slope_norm.multiply(-1).add(1)
  .focal_mean(3)
  .unitScale(0, 1);

// Lowland emphasis
var lowland = ee.Image(1).subtract(
  elev_scaled.multiply(0.6)
    .add(slope_norm.multiply(0.4))
);


// **********
//    Riparian Buffers (30m / 120m / 240m)
// **********

var buffer_30m  = flowlines.map(function(f){ return f.buffer(30);  });
var buffer_120m = flowlines.map(function(f){ return f.buffer(120); });
var buffer_240m = flowlines.map(function(f){ return f.buffer(240); });

var buffer_240m_union = buffer_240m.union();

var buffer_240m_mask = ee.Image()
  .paint(buffer_240m_union, 1)
  .toByte()
  .clip(wisc_roi);


// **********
//    Stream Proximity Weight
// **********

var dist_img = ee.Image().toFloat()
  .paint(flowlines, 1)
  .fastDistanceTransform(30)
  .sqrt()
  .clip(wisc_roi);

var prox_weight = dist_img.expression(
  'exp(-b * (d / 30))',
  {d: dist_img, b: 0.2}
);


// **********
//    Final NRI based on year of interest.
// **********

var NRI_2023 = ndvi_scaled
  .multiply(prox_weight)
  .multiply(lowland)
  .multiply(flow_proxy)
  .updateMask(buffer_240m_mask)
  .clip(wisc_roi)
  .unitScale(0, 1)
  .rename('NRI_Wisconsin');


// ********** Map + Export

Map.addLayer(
  NRI_2023,
  {min: 0, max: 1,
   palette: ['white','#e0f3db','#a8ddb5','#43a2ca','#0868ac']},
  'NRI 2023 — Wisconsin'
);

Map.addLayer(
  ee.Image().paint(allWatersheds, 1, 2),
  {palette: ['black']},
  'HUC12 Watersheds',
  false
);

Export.image.toDrive({
  image: NRI_2023,
  description: 'NRI_2023_Wisconsin_RiparianMasked',
  folder: 'NRI_Final',
  fileNamePrefix: 'NRI_2023_Wisconsin_RiparianMasked',
  region: wisc_roi.geometry(),
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Map.centerObject(wisc_roi, 9);

// ***************************************************************
// Partner Module — Watershed => Water Quality => TNRI     
//
// Instructions:
//
// 1. Upload your water-quality dataset as a TABLE asset:
//       Code Editor => Assets => NEW => Table upload
//
//    Required columns:
//       1 'system:time_start'   (timestamp; can be monthly, weekly, etc.)
//       2 'value'               (nutrient concentration, mg/L; dissolved oxygen, mg/L/saturation; turbidity, FNU; any indicator parameter works)
//
//      example CSV:
//        system:time_start,value
//        2023-03-01T00:00:00Z,1.84
//        2023-04-15T00:00:00Z,2.10
//        2023-05-01T00:00:00Z,1.95
//
//    After upload, copy your asset ID and replace it below.
//
// 2. Click inside any watershed on the map.
//
//    The script will:
//       - Identify watershed boundaries
//       - Clip the NRI surface to that watershed
//       - Compute the mean nutrient concentration
//       - Multiply by YOUR specified discharge
//       - Export TNRI (a GeoTIFF) to Google Drive


// ********** 1. Identify watershed from map click **********
function getWatershedFromClick(pt) {
  return allWatersheds.filterBounds(pt).first();
}

Map.onClick(function(coords) {

  var pt = ee.Geometry.Point(coords.lon, coords.lat);
  print('Selected point:', pt);

  var ws = getWatershedFromClick(pt);
  print('Identified Watershed:', ws);

  var ws_geom = ee.Feature(ws).geometry();


  // ********** 3. Clip base NRI to the selected watershed 
  var wsNRI = NRI_2023.clip(ws_geom);

  Map.addLayer(
    wsNRI,
    {min: 0, max: 1,
     palette: ['white','#e0f3db','#a8ddb5','#43a2ca','#0868ac']},
    'NRI for Selected Watershed'
  );


  // ********** 4. Load Partner Water Quality (Add data according to your .csv file nomenclature)
  //
  // Replace this path with your water quality feature collection.
  // Column "value" must contain nutrient concentration (mg/L).
  //
  var partner_WQ = ee.FeatureCollection(
    'users/partner_username/WaterQuality_Data'
  );


  // ********** 5. Compute MEAN nutrient concentration (mg/L) 
  // This uses simple averaging logic - friendly to continuous or non-continuous in situ data.
  
  var mean_WQ = partner_WQ.aggregate_mean('value');


  // ********** 6. Partner select discharge according to USGS stream gauge & format from cfs to (m^3/s) 
  //
  // Instructions:
  // - Go to: https://waterdata.usgs.gov/nwis/uv (Or search for nearest USGS stream gauge manually e.g. https://waterdata.usgs.gov/monitoring-location/USGS-05545750/#dataTypeId=continuous-00060-0&period=P7D&showFieldMeasurements=true)
  // - Select your state, discharge (cfs) nearest stream gauge
  // - Obtain mean discharge for the same time window as your WQ data
  // - Convert from cubic feet per second (cfs) to m^3/s:
  //
  //       m^3/s = cfs * 0.0283168
  //
  // - Enter that value below by replacing '1'.
  //
  var discharge = ee.Number(1);   // <── EDIT THIS (1) VALUE ONLY


  // ********** 7. Compute nutrient flux (mg/s) 
  //
  // (mg/L) * (m^3/s) * (1000) (L/m^3) => mg/s   (* = multiplication)
  //
  var flux = mean_WQ.multiply(discharge).multiply(1000);
  var flux_img = ee.Image.constant(flux).rename('flux');


  // ********** 8. MERIT Hydro for localized hydrology conditions 
  //
  var merit = ee.Image('MERIT/Hydro/v1_0_1')
    .select('upa')
    .clip(ws_geom);

  var flow_acc_norm = merit
    .log10()
    .unitScale(0, 3)
    .pow(2)
    .clamp(0, 1)
    .rename('flow_acc_norm');


// ********** 9. Compute TNRI 
//
// Formula:
//
//   TNRI = (NRI) * (MERIT_hydro) * (nutrient_flux)
//
  var TNRI = wsNRI
    .multiply(flow_acc_norm)
    .multiply(flux_img)
    .rename('TNRI');

   Map.addLayer(
     TNRI,
     {min: 0, max: 5e6,
      palette: ['white','lightgreen','green','darkgreen']},
     'TNRI Selected'
   );


  // ********** 10. Export TNRI to Google Drive **********
  Export.image.toDrive({
    image: TNRI,
    description: 'TNRI_Selected_Watershed',
    region: ws_geom,
    scale: 30,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });

});

Hydro_Map_Base.js
/Users/aaronwilson/Downloads/Hydro_Map_Base.js
