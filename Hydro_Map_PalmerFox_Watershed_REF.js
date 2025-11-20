// *************** ROI + Tables 
var county_list = ["Door","Kewaunee","Manitowoc","Sheboygan","Ozaukee","Milwaukee","Racine","Kenosha"];
var wisc_roi = table.filter(ee.Filter.inList("COUNTY_NAM", county_list));
var palmerfox = table2.filter(ee.Filter.eq('HUC12_CODE', '071200061003'));
Map.addLayer(ee.Image().byte().paint({featureCollection: palmerfox, color: 1, width: 2}),
             {palette: 'red'}, 'Palmer Creek–Fox River Watershed Boundary');
var wdnr_flowlines_clipped = table3.filterBounds(palmerfox);

// *************** CDL Winter Cover Crops of Interest
var winter_crops = ee.List([
  24,26,27,36,224,225,226,228,230,231,232,233,234,235,236,237,238,239,240,241,254
]).distinct();

var winter_vals = ee.List.repeat(1, winter_crops.length());

// *************** CDL Collection for all available years 
var cdl_coll = ee.ImageCollection("USDA/NASS/CDL")
  .filterBounds(wisc_roi)
  .filterDate("2008-01-01", "2024-12-31")
  .select("cropland")
  .sort("system:time_start");

// *************** Landsat helpers
function maskLandsatL2(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0) // clouds
    .and(qa.bitwiseAnd(1 << 1).eq(0))    // dilated
    .and(qa.bitwiseAnd(1 << 4).eq(0))    // shadows
    .and(qa.bitwiseAnd(1 << 5).eq(0));   // snow
  return img.updateMask(mask);
}
function scaleSR(img) {
  var srBands = img.bandNames().filter(ee.Filter.stringContains('item', 'SR_B'));
  var scaled  = img.select(srBands).multiply(0.0000275).add(-0.2);
  return img.addBands(scaled, null, true);
}
function addNDVI(img) {
  var bands = img.bandNames();
  var hasB5 = bands.contains('SR_B5');
  var nir = ee.String(ee.Algorithms.If(hasB5, 'SR_B5', 'SR_B4'));
  var red = ee.String(ee.Algorithms.If(hasB5, 'SR_B4', 'SR_B3'));
  var ndvi = img.normalizedDifference([nir, red]).rename('NDVI');
  return img.addBands(ndvi);
}
function getSpringNDVI(year) {
  var start = ee.Date.fromYMD(year, 2, 1);
  var end   = ee.Date.fromYMD(year, 5, 31);
  var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2");
  var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2");
  var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
  var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");
  var coll = ee.ImageCollection(
    ee.Algorithms.If(year <= 2011, l5,
      ee.Algorithms.If(year === 2012, l7, l8.merge(l9)))
  )
    .filterBounds(wisc_roi)
    .filterDate(start, end)
    .map(maskLandsatL2)
    .map(scaleSR)
    .map(addNDVI)
    .select('NDVI');
  return coll.median().clip(wisc_roi);
}

// *************** Spring NDVI 2023 
var springNDVI_2023 = getSpringNDVI(2023);
var springNDVI_2023_clipped = springNDVI_2023.clip(palmerfox.geometry());
Map.addLayer(
  springNDVI_2023_clipped,
  {min: 0, max: 0.8, palette: ['#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#31a354','#006837']},
  'Spring NDVI 2023 (Clipped)', true
);
Export.image.toDrive({
  image: springNDVI_2023_clipped,
  description: 'Spring_NDVI_2023_Clipped',
  folder: 'NRI_Final',
  fileNamePrefix: 'Spring_NDVI_2023_Clipped',
  region: palmerfox.geometry(),
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// *************** USGS Stream Gauge point
// Replicate line 87 with your own coordinates
// of where in-situ water quality data was collected 
var gauge_point = ee.Geometry.Point([-88.225923, 42.610851]);
Map.addLayer(gauge_point, {color: 'red'}, 'Fox River Stream Gauge');

// *************** RIPARIAN BUFFERS (30/120/240) + visualize
var buffer_30m  = wdnr_flowlines_clipped.map(function(f){ return f.buffer(30);  });
var buffer_120m = wdnr_flowlines_clipped.map(function(f){ return f.buffer(120); });
var buffer_240m = wdnr_flowlines_clipped.map(function(f){ return f.buffer(240); });
var buffer_30m_union  = buffer_30m.union();
var buffer_120m_union = buffer_120m.union();
var buffer_240m_union = buffer_240m.union();
Map.addLayer(wdnr_flowlines_clipped, {color: 'red'},      'WI DNR Flowlines (Clipped)');
Map.addLayer(buffer_30m_union,       {color: 'lightblue'}, '30m Riparian Buffer (Clipped)');
Map.addLayer(buffer_120m_union,      {color: 'blue'},      '120m Riparian Buffer (Clipped)');
Map.addLayer(buffer_240m_union,      {color: 'darkblue'},  '240m Riparian Buffer (Clipped)');

// *************** Proximity weight (exponential decay from streams)
var dist_img = ee.Image().toFloat()
  .paint(wdnr_flowlines_clipped, 1)
  .fastDistanceTransform(30).sqrt()
  .clip(palmerfox.geometry());
var prox_weight = dist_img.expression('exp(-b * d)', {
  d: dist_img.divide(30),  // pixels
  b: 0.2
}).rename('prox_w');
var prox_norm = prox_weight.updateMask(prox_weight.gte(0.01));

// *************** SRTM + Slope normalizations (for hydro/terrain contexts)
var dem   = ee.Image("USGS/SRTMGL1_003").clip(palmerfox.geometry());
var slope = ee.Terrain.slope(dem);
var ndvi_scaled_2023 = springNDVI_2023.unitScale(0, 0.8);
var slope_norm = slope.unitScale(0, 15);
var slope_flat = slope_norm.expression('1 - s', {s: slope_norm});

// *************** NRI 
var elev_scaled  = dem.unitScale(100, 400).clamp(0, 1);
var flow_proxy   = slope_norm.multiply(-1).add(1).focal_mean(3).unitScale(0, 1);
var lowland_weight = ee.Image(1).subtract(elev_scaled.multiply(0.6).add(slope_norm.multiply(0.4)));
var hydro_weight = lowland_weight.multiply(flow_proxy);
var riparian_mask_geom = wdnr_flowlines_clipped.map(function(f){ return f.buffer(240); }).union();
var riparian_mask_img  = ee.Image().paint(riparian_mask_geom, 1).toByte();
var weight = hydro_weight.multiply(prox_weight);
var NRI_2023_SRTM = ndvi_scaled_2023.multiply(weight)
  .updateMask(riparian_mask_img)
  .clip(palmerfox.geometry())
  .rename('NRI_2023_SRTM');
Map.addLayer(
  NRI_2023_SRTM,
  {min: 0, max: 1, palette: ['white','#e0f3db','#a8ddb5','#43a2ca','#0868ac']},
  'NRI (SRTM 30m)', true
);
Export.image.toDrive({
  image: NRI_2023_SRTM,
  description: 'NRI_2023_PalmerFox',
  folder: 'NRI_Final',
  fileNamePrefix: 'NRI_2023_PalmerFox',
  region: palmerfox.geometry(),
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// *************** MERIT HYDRO 
var flow_acc = ee.Image('MERIT/Hydro/v1_0_1').select('upa').clip(palmerfox.geometry());
var flow_acc_norm = flow_acc.log10().unitScale(0, 3).pow(2).clamp(0, 1).rename('flow_acc_norm');

// *************** Nitrate + Nitrite example input
var nitrateNitrite  = ee.FeatureCollection('projects/fbrtransport/assets/Nitrate_Nitrite_Fox'); //Import WQ data here
function toThreeHourMean(fc) {
  var withBucket = fc.map(function(f) {
    var t = ee.Date.parse("YYYY-MM-dd'T'HH:mm:ss'Z'", f.getString('system:time_start'));
    var bucketMillis = t.millis().divide(1000 * 60 * 60 * 3).floor().multiply(1000 * 60 * 60 * 3);
    return f.set({bucketMillis: bucketMillis, value: ee.Number.parse(f.get('value'))});
  });
  var groups = ee.List(
    withBucket.reduceColumns({
      selectors: ['bucketMillis', 'value'],
      reducer: ee.Reducer.mean().group({groupField: 0, groupName: 'bucketMillis'})
    }).get('groups')
  );
  return ee.FeatureCollection(groups.map(function(g) {
    g = ee.Dictionary(g);
    var bm = ee.Number(g.get('bucketMillis'));
    var meanVal = ee.Number(g.get('mean'));
    return ee.Feature(null, {date: ee.Date(bm), value: meanVal});
  })).sort('date');
}
var nitrate3hr = toThreeHourMean(nitrateNitrite);
var nitrateChart = ui.Chart.feature.byFeature({
  features: nitrate3hr, xProperty: 'date', yProperties: ['value']
}).setChartType('LineChart').setOptions({
  title: 'Fox River Nitrate + Nitrite (2023)', lineWidth: 1, pointSize: 0, dataOpacity: 1.0,
  colors: ['blue'], hAxis: {title: 'Date', format: 'MMM d, yyyy', gridlines: {count: 8}},
  vAxis: {title: 'Value (mg/L as N)'}, legend: {position: 'none'}
});
print(nitrateChart);
// Edit start and end date to reflect study period of interest 
var startDate = '2023-02-01';
var endDate   = '2023-05-31';
var table6 = nitrateNitrite.map(function(f){  // ensure time_start recognized
  return f.set('system:time_start', ee.Date(f.get('system:time_start')));
});
var nitrate_filtered = table6.filter(ee.Filter.date(startDate, endDate));
var mean_nitrate   = nitrate_filtered.aggregate_mean('value');
print('Mean Nitrate (mg/L):', mean_nitrate);

// *************** Constants for this example 
var meanNitrate   = 2.13;      // mg/L (excel)
var meanDischarge = 1145.68;   // m^3/s (gage). This is important to collect from USGS stream gauge online data export

// *************** Spatialized nitrate load (base + MERIT-weighted)
var nitrateFlux_const = ee.Image.constant(meanNitrate * meanDischarge * 1000).rename('nitrate_flux_mgs');
var ndvi_norm = ndvi_scaled_2023.updateMask(ndvi_scaled_2023.gte(0.2));
var nitrateRetentionProxy = ndvi_norm
  .multiply(prox_norm)
  .multiply(slope_flat)
  .rename('nitrate_retention_proxy');
var spatialNitrateLoad = nitrateFlux_const.multiply(nitrateRetentionProxy)
  .rename('nitrate_load_spatialized')
  .clip(buffer_240m_union);
Map.addLayer(spatialNitrateLoad,
  {min: 1, max: 3e6, palette: ['white', 'lightgreen', 'green']},
  'Spatialized Nitrate Load (mg/s)'
);

// MERIT-weighted (stronger hydrology influence)
var nitrateFlux_merit = flow_acc_norm.multiply(meanNitrate * meanDischarge * 1000)
  .rename('nitrate_flux_mgs');
var ndvi_weighted = ndvi_scaled_2023.pow(0.5);
var nitrateRetentionProxy_MERIT = ndvi_weighted
  .multiply(prox_weight)
  .multiply(slope_flat)
  .multiply(flow_acc_norm)
  .rename('nitrate_retention_proxy_merit');
var spatialNitrateLoad_MERIT = nitrateFlux_merit.multiply(nitrateRetentionProxy_MERIT)
  .rename('nitrate_load_spatialized_merit')
  .clip(palmerfox.geometry());
Map.addLayer(spatialNitrateLoad_MERIT,
  {min: 1, max: 3e6, palette: ['white', '#c7e9b4', '#7fcdbb', '#1d91c0', '#0c2c84']},
  'Spatialized Nitrate Load (Improved MERIT flow weighting, mg/s)'
);

// *************** Watershed-wide leaching and buildup 
var moderate_flow = flow_acc_norm.gt(0.2).and(flow_acc_norm.lt(0.6));
var near_stream   = prox_norm.gte(0.3);
var nitrate_leach_potential = ndvi_norm
  .multiply(slope_flat)
  .multiply(moderate_flow)
  .multiply(near_stream)
  .multiply(nitrateFlux_merit.divide(1e6))
  .rename('nitrate_leach_potential');
var nitrate_leach_norm = nitrate_leach_potential.unitScale(0, 1).clamp(0, 1)
  .clip(palmerfox.geometry());
Map.addLayer(nitrate_leach_norm, {
  min: 0, max: 1, palette: ['#ffffcc','#c2e699','#78c679','#31a354','#006837']
}, 'Nitrate Leaching Potential (Watershed-wide)');

var hydro_context = flow_acc_norm.pow(0.25);
var proximity_decay_build = prox_weight.expression('exp(-b*d)', {d: prox_weight.divide(60), b: 0.2});
var nitrate_buildup_potential = ndvi_norm
  .multiply(slope_flat)
  .multiply(proximity_decay_build)
  .multiply(hydro_context)
  .rename('nitrate_buildup_potential');
var nitrate_buildup_clipped = nitrate_buildup_potential.unitScale(0, 1)
  .clip(palmerfox.geometry());
Map.addLayer(nitrate_buildup_clipped, {
  min: 0, max: 1, palette: ['#ffffe5','#d9f0a3','#78c679','#238443','#004529']
}, 'Nitrate Buildup Potential Watershed-wide');

// *************** Terrestrial nitrate buildup (main)
var low_flow_focus = flow_acc_norm.expression('1 - f', {f: flow_acc_norm}).clamp(0, 1);
var proximity_decay_land = prox_weight.expression('exp(-b*d)', {d: prox_weight.divide(90), b: 0.15});
var nitrate_land_buildup = ndvi_norm
  .multiply(slope_flat)
  .multiply(low_flow_focus)
  .multiply(proximity_decay_land)
  .rename('nitrate_land_buildup');
var nitrate_land_buildup_adj = nitrate_land_buildup
  .pow(1.5)
  .multiply(slope_norm.expression('1 - (s * 1.2)', {s: slope_norm}))
  .clamp(0, 1);
Map.addLayer(nitrate_land_buildup_adj, {
  min: 0, max: 0.7,
  palette: ['white','#e0f3db','#a8ddb5','#43a2ca','#0868ac']
}, 'NitrogenCover_2023');
Export.image.toDrive({
  image: nitrate_land_buildup_adj,
  description: 'NitrogenCover_2023',
  folder: 'NRI_Final',
  fileNamePrefix: 'NitrogenCover_2023',
  region: palmerfox.geometry(),
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// *************** 2023 CDL winter crops example. Change year as needed 
var cdl_2023 = ee.Image(cdl_coll.filter(ee.Filter.calendarRange(2023, 2023, "year")).first());
var winter_mask = cdl_2023.remap({
  from: winter_crops, to: winter_vals, defaultValue: 0, bandName: "cropland"
}).rename("winter_mask").clip(wisc_roi);

// *************** Terrestrial nitrogen retention cover areas highlighting overlap with CDL winter crop presence
var TNRP = nitrate_land_buildup_adj.select(0).rename('TNRP').clip(palmerfox.geometry());
var cover_mask_vis = winter_mask.clip(palmerfox.geometry());
var TNRP_cover    = TNRP.updateMask(cover_mask_vis.eq(1));

Map.addLayer(TNRP_cover,    {min: 0, max: 0.3, palette: ['lightyellow','orange','red']}, 'TNRP (Cover Areas)');


Export.image.toDrive({
  image: TNRP_cover,
  description: 'TNRP_CoverAreas_2023',
  folder: 'NRI_Final',
  fileNamePrefix: 'TNRP_CoverAreas_2023',
  region: palmerfox.geometry(),
  scale: 30,
  maxPixels: 1e13
});

// *************** Buffer stats (30–500 m) for data analysis
var buffer_360m = wdnr_flowlines_clipped.map(function(f){ return f.buffer(360); });
var buffer_500m = wdnr_flowlines_clipped.map(function(f){ return f.buffer(500); });
var buffer_360m_union = buffer_360m.union();
var buffer_500m_union = buffer_500m.union();

function tagBuffer(fc, dist) { return fc.map(function(f){ return f.set('Buffer_m', dist); }); }
var buffers_all = ee.FeatureCollection([
  tagBuffer(buffer_30m, 30),
  tagBuffer(buffer_120m, 120),
  tagBuffer(buffer_240m, 240),
  tagBuffer(buffer_360m, 360),
  tagBuffer(buffer_500m, 500)
]).flatten();

// combine NDVI, NRI, and MERIT load for stats
var combo = springNDVI_2023.select([0]).rename('NDVI')
  .addBands(NRI_2023_SRTM.rename('NRI'))
  .addBands(spatialNitrateLoad_MERIT.rename('Nitrate_Load'));

var stats = combo.reduceRegions({
  collection: buffers_all,
  reducer: ee.Reducer.mean(),
  scale: 30
});
Export.table.toDrive({
  collection: stats,
  description: 'BufferStats_NDVI_NRI_NitrateLoad_30to500m',
  folder: 'NRI_Final',
  fileFormat: 'CSV'
});

// *************** View
Map.centerObject(palmerfox, 12);

Hydro_Map_PalmerFox_Watershed_REF.js
/Users/aaronwilson/Downloads/Hydro_Map_PalmerFox_Watershed_REF.js