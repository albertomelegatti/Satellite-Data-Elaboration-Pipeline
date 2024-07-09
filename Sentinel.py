from sentinelhub import SentinelHubRequest, DataCollection, MimeType, CRS, BBox, SHConfig, Geometry
import sys
import pyproj
import os
from function import getValuesMatrix, getPercentageMatrix, generateImage, compareDates, compareCoordinates, extractImagesFromTar, jsonBuilder, convertImageTo4Colors, htmlLeafLetBuilder, tiff2Jpg
import datetime

# Verifica che il numero di argomenti sia corretto
if sys.argv[10] == "gprox":
    if len(sys.argv) != 12 and len(sys.argv) != 13:
        print("Usage: python mainScript.py <client_id> <client_secret> <NW Latitude> <NW Longitude> <SE Latitude> <Se Longitude> <timeIntervalStart> <timeIntervalEnd> <maxCloudCoverage> <product> <metersRadius> (<html>)")
        sys.exit(1)

else:
    if len(sys.argv) != 11 and len(sys.argv) != 12:
        print("Usage: python mainScript.py <client_id> <client_secret> <NW Latitude> <NW Longitude> <SE Latitude> <SE Longitude> <timeIntervalStart> <timeIntervalEnd> <maxCloudCoverage> <product> (<html>)")
        sys.exit(1)


# Verifica che la data di inizio sia precedente alla data di fine
if compareDates(sys.argv[7], sys.argv[8]) == 2:
    print("The start date must be before the end date")
    sys.exit(1)

# Verifica che le coordinate siano corrette
if compareCoordinates(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])) == 1:
    print("The coordinates must be within the following ranges: latitude [-90, 90], longitude [-180, 180]")
    sys.exit(1)
elif compareCoordinates(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])) == 2:
    print("The coordinates must be entered in the correct order: NW_Latitude > SE_Latitude, NW_Longitude < SE_Longitude")
    sys.exit(1)

# Verifica che il prodotto scelto sia tra quelli disponibili
if sys.argv[10] != "stemp" and sys.argv[10] != "stype" and sys.argv[10] != "gprox" and sys.argv[10] != "vis":
    print("Invalid choice")
    sys.exit(1)

# Verifica che il valore di maxCloudCoverage sia compreso tra 0 e 100
if float(sys.argv[9]) < 0 or float(sys.argv[9]) > 100:
    print("The value of maxCloudCoverage must be between 0 and 100")
    sys.exit(1)


htmlFlag = "false"

# Nel caso in cui il prodotto sia gprox, viene richiesto un parametro in più, metersRadius
if sys.argv[10] == "gprox":
    
    # Assegnamento dei parametri alle variabili
    client_id = sys.argv[1]
    client_secret = sys.argv[2]
    NW_Latitude = float(sys.argv[3])
    NW_Longitude = float(sys.argv[4])
    SE_Latitude = float(sys.argv[5])
    SE_Longitude = float(sys.argv[6])
    timeIntervalStart = sys.argv[7]
    timeIntervalEnd = sys.argv[8]
    maxCloudCoverage = float(sys.argv[9])
    product = "gprox"
    metersRadius = int(sys.argv[11])
    
    if len(sys.argv) == 13 and sys.argv[12] == "html":
        htmlFlag = "true"
    
    # Verifica che il raggio sia maggiore di 0
    if metersRadius <= 0:
        print("Radius must be greater than 0")
        sys.exit(1)


else:
    # Assegnamento dei parametri alle variabili
    client_id = sys.argv[1]
    client_secret = sys.argv[2]
    NW_Latitude = float(sys.argv[3])
    NW_Longitude = float(sys.argv[4])
    SE_Latitude = float(sys.argv[5])
    SE_Longitude = float(sys.argv[6])
    timeIntervalStart = sys.argv[7]
    timeIntervalEnd = sys.argv[8]
    maxCloudCoverage = float(sys.argv[9])
    product = sys.argv[10]
    
    
    if len(sys.argv) == 12 and sys.argv[11] == "html":
        htmlFlag = "true"



# Calcolo della risoluzione
# Vengono calcolate le dimensioni in metri dell'area di interesse
# Vengono calcolate le dimensioni in pixel dell'area di interesse
# Viene calcolata la risoluzione (dimensioni in metri / dimensioni in pixel)
wgs84 = pyproj.Geod(ellps='WGS84')

lat1, lon1 = NW_Latitude, NW_Longitude
lat2, lon2 = SE_Latitude, SE_Longitude

verticalSideMeter = wgs84.inv(lon1, lat1, lon1, lat2)[2]
horizontalSideMeter = wgs84.inv(lon1, lat1, lon2, lat1)[2]

longSide = max(verticalSideMeter, horizontalSideMeter)
shortSide = min(verticalSideMeter, horizontalSideMeter)

verticalSidePixel = verticalSideMeter * (750 / longSide)
horizontalSidePixel = horizontalSideMeter * (750 / longSide)

resolution = longSide / 750

if product == "gprox":
    # Il valore di metri di raggio sono i metri di raggio che l'utente vuole considerare per il calcolo dell'indice di prossimità di verde
    # Dividendo i metri di raggio per la risoluzione si ottiene il raggio in pixel, che verrà utilizzato per il calcolo dell'indice di prossimità di verde
    vegetationIndexRadius = (float(metersRadius) / resolution)
    vegetationIndexRadius = (metersRadius / resolution)

    # Arrotondamento del raggio per l'indice di prossimità di verde
    vegetationIndexRadius = round(vegetationIndexRadius)


# Creazione della cartella "results" in cui verranno salvati tutti gli output
folderPath = os.path.dirname(os.path.realpath(__file__))
resultFolderName = "results_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") 
resultFolderPath = os.path.join(folderPath, resultFolderName)

    
# Controlla se la cartella esiste già
if not os.path.exists(resultFolderPath):
    os.makedirs(resultFolderPath)

percentage_matrix = []


# Determina quale script eseguire in base alla scelta
# LAND SURFACE TEMPERATURE
if product == 'stemp':
    
    config = SHConfig()
    config.sh_client_id = client_id
    config.sh_client_secret = client_secret
    evalscript = """
       // VERSION 3

/**
  This script is directly based on the Landsat-8 Land Surface Temperature Mapping script by Mohor Gartner
  https://custom-scripts.sentinel-hub.com/landsat-8/land_surface_temperature_mapping/
  since the script uses Landsat TIRS B10 for brightness temperature 
  mapping and Landsat OLI NDVI to scale for emissivity, this can be followed using 
  Sentinel-3 SLSTR S08 and Sentinel-3 OLCI NDVI

  in order to use this script you have to enable "use additional datasets (advanced)"
  and set S-3 OLCI and S-3 SLSTR as the primary and additional dataset.

  Aliases should be 
   - Sentinel-3 OLCI=S3OLCI
   - Sentinel-3 SLSTR=S3SLSTR

  STARTING OPTIONS
  for analysis of one image (EO Browser), choose option=0. In case of MULTI-TEMPORAL analyis, 
  option values are following:
  0 - outputs average LST in selected timeline (% of cloud coverage should be low, e.g. < 10%)
  1 - outputs maximum LST in selected timeline (% of cloud coverage can be high)
  2 - THIS OPTION IS CURRENTLY NOT FUNCTIONAL - outputs standard deviation LST in selected timeline; 
      minTemp and highTemp are overwritten with values 0 and 10 (% of cloud coverage should be low, e.g. < 5%)
*/
var option = 0;

// minimum and maximum values for output colour chart red to white for temperature in °C. Option 2 overwrites this selection!
var minC = 0;
var maxC = 50;


////INPUT DATA - FOR BETTER RESULTS, THE DATA SHOULD BE ADJUSTED
// NVDIs for bare soil and NDVIv for full vegetation
// Note: NVDIs for bare soil and NDVIv for full vegetation are needed to 
//       be evaluated for every scene. However in the custom script, default values are set regarding:
// https://profhorn.meteor.wisc.edu/wxwise/satmet/lesson3/ndvi.html 
// https://www.researchgate.net/post/Can_anyone_help_me_to_define_a_range_of_NDVI_value_to_extract_bare_soil_pixels_for_Landsat_TM
// NVDIs=0.2, NDVIv=0.8
// other source suggests global values: NVDIs=0.2, NDVIv=0.5; 
// https://www.researchgate.net/publication/296414003_Algorithm_for_Automated_Mapping_of_Land_Surface_Temperature_Using_LANDSAT_8_Satellite_Data
var NDVIs = 0.2;
var NDVIv = 0.8;

// emissivity
var waterE = 0.991;
var soilE = 0.966;
var vegetationE = 0.973;
//var buildingE=0.962;
var C = 0.009; //surface roughness, https://www.researchgate.net/publication/331047755_Land_Surface_Temperature_Retrieval_from_LANDSAT-8_Thermal_Infrared_Sensor_Data_and_Validation_with_Infrared_Thermometer_Camera

//central/mean wavelength in meters, Sentinel-3 SLSTR B08 (almost the same as Landsat B10)
var bCent = 0.000010854;

// rho =h*c/sigma=PlanckC*velocityLight/BoltzmannC
var rho = 0.01438; // m K

//// visualization
// if result should be std dev (option=2), overwrite minMaxC.
if (option == 2) {
  minC = 0;
  maxC = 25;
}
let viz = ColorGradientVisualizer.createRedTemperature(minC, maxC);

//this is where you set up the evalscript to access the bands of the two datasets in the fusion

function setup() {
  return {
    input: [
      { datasource: "S3SLSTR", bands: ["S8"] },
      { datasource: "S3OLCI", bands: ["B06", "B08", "B17"] }],
    output: [
      { id: "default", bands: 3, sampleType: SampleType.AUTO }
    ],
    mosaicking: "ORBIT"
  }
}


//emissivity calc (Unchanged from Landsat script)
//https://www.researchgate.net/publication/296414003_Algorithm_for_Automated_Mapping_of_Land_Surface_Temperature_Using_LANDSAT_8_Satellite_Data
//https://www.academia.edu/27239873/Investigating_Land_Surface_Temperature_Changes_Using_Landsat_Data_in_Konya_Turkey
function LSEcalc(NDVI, Pv) {
  var LSE;
  if (NDVI < 0) {
    //water
    LSE = waterE;
  } else if (NDVI < NDVIs) {
    //soil
    LSE = soilE;
  } else if (NDVI > NDVIv) {
    //vegetation
    LSE = vegetationE;
  } else {
    //mixtures of vegetation and soil
    LSE = vegetationE * Pv + soilE * (1 - Pv) + C;
  }
  return LSE;
}

function evaluatePixel(samples) {
  // starting values max, avg, stdev, reduce N, N for multi-temporal
  var LSTmax = -999;
  var LSTavg = 0;
  var LSTstd = 0;
  var reduceNavg = 0;
  var N = samples.S3SLSTR.length;

  //to caputure all values of one pixel for for whole timeline in mosaic order
  var LSTarray = [];

  // multi-temporal: loop all samples in selected timeline
  for (let i = 0; i < N; i++) {
    //// for LST S8
    var Bi = samples.S3SLSTR[i].S8;
    var B06i = samples.S3OLCI[i].B06;
    var B08i = samples.S3OLCI[i].B08;
    var B17i = samples.S3OLCI[i].B17;

    // some images have errors, whole area is either B10<173K or B10>65000K. Also errors, where B06 and B17 =0. Therefore no processing if that happens, in addition for average and stdev calc, N has to be reduced!
    if ((Bi > 173 && Bi < 65000) && (B06i > 0 && B08i > 0 && B17i > 0)) {
      // ok image
      //1 Kelvin to C
      var S8BTi = Bi - 273.15;
      //2 NDVI - Normalized Difference vegetation Index - based on this custom script: https://custom-scripts.sentinel-hub.com/sentinel-3/ndvi/
      var NDVIi = (B17i - B08i) / (B17i + B08i);
      //3 PV - proportional vegetation
      var PVi = Math.pow(((NDVIi - NDVIs) / (NDVIv - NDVIs)), 2);
      //4 LSE land surface emmisivity  
      var LSEi = LSEcalc(NDVIi, PVi);
      //5 LST
      var LSTi = (S8BTi / (1 + (((bCent * S8BTi) / rho) * Math.log(LSEi))));

      ////temporary calculation
      //avg
      LSTavg = LSTavg + LSTi;
      //max
      if (LSTi > LSTmax) { LSTmax = LSTi; }
      //array
      LSTarray.push(LSTi);
    } else {
      // image NOT ok
      ++reduceNavg;
    }
  }
  // correct N value if some images have errors and are not analysed
  N = N - reduceNavg;

  // calc final avg value
  LSTavg = LSTavg / N;

  // calc final stdev value
  for (let i = 0; i < LSTarray.length; i++) {
    LSTstd = LSTstd + (Math.pow(LSTarray[i] - LSTavg, 2));
  }
  LSTstd = (Math.pow(LSTstd / (LSTarray.length - 1), 0.5));

  // WHICH LST to output, it depends on option variable: 0 for one image analysis (OE Browser); MULTI-TEMPORAL: 0->avg; 1->max; 2->stdev
  let outLST = (option == 0)
    ? LSTavg
    : (option == 1)
      ? LSTmax
      : LSTstd;

  //// output to image
  return viz.process(outLST);
}
    """

    bbox = BBox(bbox=[NW_Longitude, SE_Latitude, SE_Longitude, NW_Latitude], crs=CRS.WGS84)
    geometry = Geometry(
        geometry={
            "coordinates": [[[NW_Longitude, SE_Latitude], [NW_Longitude, NW_Latitude], [SE_Longitude, NW_Latitude], [SE_Longitude, SE_Latitude], [NW_Longitude, SE_Latitude]]],
            "type": "Polygon"
        },
        crs=CRS.WGS84
    )
    request = SentinelHubRequest(
        evalscript = evalscript,
        data_folder = resultFolderPath,  
        input_data=[
            SentinelHubRequest.input_data(
                data_collection=DataCollection.SENTINEL3_SLSTR,
                identifier="S3SLSTR",    
                time_interval=(timeIntervalStart, timeIntervalEnd),
                other_args={"dataFilter": {"maxCloudCoverage": maxCloudCoverage}}
            ),
            SentinelHubRequest.input_data(
                data_collection=DataCollection.SENTINEL3_OLCI,
                identifier="S3OLCI",   
                time_interval=(timeIntervalStart, timeIntervalEnd),
                other_args={"dataFilter": {"maxCloudCoverage": maxCloudCoverage}}
            ),
        ],
        responses=[
            SentinelHubRequest.output_response('default', MimeType.TIFF),
            SentinelHubRequest.output_response('default', MimeType.JPG),
        ],
        bbox=bbox,
        geometry=geometry,
        size=[horizontalSidePixel, verticalSidePixel],
        config=config
    )
    print("Request sent")
    response = request.get_data(save_data=True)
    print("Request completed")

    # Estrazione delle immagini dal file tar
    extractImagesFromTar(resultFolderPath)
    


# SOIL TYPE
elif product == 'stype':
    config = SHConfig()
    config.sh_client_id = client_id
    config.sh_client_secret = client_secret
    evalscript = """
        //VERSION=3
        function evaluatePixel(samples) {
            var NDWI = index(samples.B03, samples.B08); 
            var NDVI = index(samples.B08, samples.B04);
            var BareSoil = 2.5 * ((samples.B11 + samples.B04) - (samples.B08 + samples.B02)) / ((samples.B11 + samples.B04) + (samples.B08 + samples.B02));

            if (NDWI > 0.2) {
                return [0, 0.5, 1, samples.dataMask];
            }
            if ((samples.B11 > 0.8) || (NDVI < 0.1)) {
                return [1, 1, 1, samples.dataMask];
            }
            if (NDVI > 0.2) {
                return [0, 1, 0, samples.dataMask];
            } else {
                return [1, 0, 0, samples.dataMask];
            }
        }
        function setup() {
            return {
                input: ["B02", "B03", "B04", "B08", "B11", "dataMask"],
                output: { bands: 4 }
            };
        }
    """

    bbox = BBox(bbox=[NW_Longitude, SE_Latitude, SE_Longitude, NW_Latitude], crs=CRS.WGS84)
    geometry = Geometry(
        geometry={
            "coordinates": [[[NW_Longitude, SE_Latitude], [NW_Longitude, NW_Latitude], [SE_Longitude, NW_Latitude], [SE_Longitude, SE_Latitude], [NW_Longitude, SE_Latitude]]],
            "type": "Polygon"
        },
        crs=CRS.WGS84
    )
    request = SentinelHubRequest(
        evalscript = evalscript,
        data_folder = resultFolderPath,  
        input_data=[
            SentinelHubRequest.input_data(
                data_collection=DataCollection.SENTINEL2_L2A,     
                time_interval=(timeIntervalStart, timeIntervalEnd),
                other_args={"dataFilter": {"maxCloudCoverage": maxCloudCoverage}}
            ),
        ],
        responses=[
            SentinelHubRequest.output_response('default', MimeType.TIFF),
            SentinelHubRequest.output_response('default', MimeType.JPG),
        ],
        bbox=bbox,
        geometry=geometry,
        size=[horizontalSidePixel, verticalSidePixel],
        config=config
    )
    print("Request sent")
    response = request.get_data(save_data=True)
    print("Request completed")

    # Estrazione delle immagini dal file tar
    extractImagesFromTar(resultFolderPath)
    
    # L'immagine jpg scaricata viene convertita in un'immagine a 4 colori
    convertImageTo4Colors(os.path.join(resultFolderPath, "extracted_contents", "default.jpg"), os.path.join(resultFolderPath, "extracted_contents", "convertedDefault.tif"))
    



# VEGETATION INDEX
elif product == 'gprox':
    config = SHConfig()
    config.sh_client_id = client_id
    config.sh_client_secret = client_secret
    evalscript = """
        //VERSION=3
        function evaluatePixel(samples) {
            var NDWI = index(samples.B03, samples.B08); 
            var NDVI = index(samples.B08, samples.B04);
            var BareSoil = 2.5 * ((samples.B11 + samples.B04) - (samples.B08 + samples.B02)) / ((samples.B11 + samples.B04) + (samples.B08 + samples.B02));

            if (NDWI > 0.2) {
                return [0, 1, 0, samples.dataMask];
            }
            if ((samples.B11 > 0.8) || (NDVI < 0.1)) {
                return [1, 1, 1, samples.dataMask];
            }
            if (NDVI > 0.2) {
                return [0, 1, 0, samples.dataMask];
            } else {
                return [1, 1, 1, samples.dataMask];
            }
        }
        function setup() {
            return {
                input: ["B02", "B03", "B04", "B08", "B11", "dataMask"],
                output: { bands: 4 }
            };
        }
    """

    bbox = BBox(bbox=[NW_Longitude, SE_Latitude, SE_Longitude, NW_Latitude], crs=CRS.WGS84)
    geometry = Geometry(
        geometry={
            "coordinates": [[[NW_Longitude, SE_Latitude], [NW_Longitude, NW_Latitude], [SE_Longitude, NW_Latitude], [SE_Longitude, SE_Latitude], [NW_Longitude, SE_Latitude]]],
            "type": "Polygon"
        },
        crs=CRS.WGS84
    )
    request = SentinelHubRequest(
        evalscript = evalscript,
        data_folder = resultFolderPath,  
        input_data=[
            SentinelHubRequest.input_data(
                data_collection=DataCollection.SENTINEL2_L2A,     
                time_interval=(timeIntervalStart, timeIntervalEnd),
                other_args={"dataFilter": {"maxCloudCoverage": maxCloudCoverage}}
            ),
        ],
        responses=[
            SentinelHubRequest.output_response('default', MimeType.TIFF),
            SentinelHubRequest.output_response('default', MimeType.JPG),
        ],
        bbox=bbox,
        geometry=geometry,
        size=[horizontalSidePixel, verticalSidePixel],
        config=config
    )
    print("Request sent")
    response = request.get_data(save_data=True)
    print("Request completed")

    # Estrazione delle immagini dal file tar
    extractImagesFromTar(resultFolderPath)

    #Esecuzione di uno script per la generazione dell'immagine con le 5 bande colorate
    input_output_paths = [
        (resultFolderPath + "/extracted_contents/default.tif", resultFolderPath + "/extracted_contents/result.tiff")
    ]

    for input_path, output_path in input_output_paths:
        value_matrix = getValuesMatrix(input_path)
        #print("Matrice dei valori creata per", input_path)
            
        percentage_matrix = getPercentageMatrix(value_matrix, vegetationIndexRadius)
        #print("Matrice delle percentuali creata per", input_path)
            
        generateImage(percentage_matrix, output_path)
        #print("Immagine generata per", output_path)
    
    os.rename(resultFolderPath + "/extracted_contents/default.jpg", resultFolderPath + "/extracted_contents/default1.jpg")
    tiff2Jpg(resultFolderPath + "/extracted_contents/result.tiff", resultFolderPath+ "/extracted_contents/default.jpg")




# TRUE COLOR
elif product == 'vis':
    config = SHConfig()
    config.sh_client_id = client_id
    config.sh_client_secret = client_secret
    evalscript = """
        //VERSION=3

        function setup() {
            return {
                input: ["B04", "B03", "B02", "dataMask"],
                output: { bands: 4 }
            };
        }

        // Contrast enhance / highlight compress


        const maxR = 3.0; // max reflectance

        const midR = 0.13;
        const sat = 1.2;
        const gamma = 1.8;

        function evaluatePixel(smp) {
            const rgbLin = satEnh(sAdj(smp.B04), sAdj(smp.B03), sAdj(smp.B02));
            return [sRGB(rgbLin[0]), sRGB(rgbLin[1]), sRGB(rgbLin[2]), smp.dataMask];
        }

        function sAdj(a) {
            return adjGamma(adj(a, midR, 1, maxR));
        }

        const gOff = 0.01;
        const gOffPow = Math.pow(gOff, gamma);
        const gOffRange = Math.pow(1 + gOff, gamma) - gOffPow;

        function adjGamma(b) {
            return (Math.pow((b + gOff), gamma) - gOffPow) / gOffRange;
        }

        // Saturation enhancement

        function satEnh(r, g, b) {
            const avgS = (r + g + b) / 3.0 * (1 - sat);
            return [clip(avgS + r * sat), clip(avgS + g * sat), clip(avgS + b * sat)];
        }

        function clip(s) {
            return s < 0 ? 0 : s > 1 ? 1 : s;
        }

        //contrast enhancement with highlight compression

        function adj(a, tx, ty, maxC) {
            var ar = clip(a / maxC, 0, 1);
            return ar * (ar * (tx / maxC + ty - 1) - ty) / (ar * (2 * tx / maxC - 1) - tx / maxC);
        }

        const sRGB = (c) => c <= 0.0031308 ? (12.92 * c) : (1.055 * Math.pow(c, 0.41666666666) - 0.055);  
    """

    bbox = BBox(bbox=[NW_Longitude, SE_Latitude, SE_Longitude, NW_Latitude], crs=CRS.WGS84)
    geometry = Geometry(
        geometry={
            "coordinates": [[[NW_Longitude, SE_Latitude], [NW_Longitude, NW_Latitude], [SE_Longitude, NW_Latitude], [SE_Longitude, SE_Latitude], [NW_Longitude, SE_Latitude]]],
            "type": "Polygon"
        },
        crs=CRS.WGS84
    )
    request = SentinelHubRequest(
        evalscript = evalscript,
        data_folder = resultFolderPath,  
        input_data=[
            SentinelHubRequest.input_data(
                data_collection=DataCollection.SENTINEL2_L1C,     
                time_interval=(timeIntervalStart, timeIntervalEnd),
                other_args={"dataFilter": {"maxCloudCoverage": maxCloudCoverage}}
            ),
        ],
        responses=[
            SentinelHubRequest.output_response('default', MimeType.TIFF),
            SentinelHubRequest.output_response('default', MimeType.JPG),
        ],
        bbox=bbox,
        geometry=geometry,
        size=[horizontalSidePixel, verticalSidePixel],
        config=config
    )
    print("Request sent")
    response = request.get_data(save_data=True)
    print("Request completed")

    # Estrazione delle immagini dal file tar
    extractImagesFromTar(resultFolderPath)

else:
    print("Invalid choice")
    sys.exit(1)

# In tutti i casi tranne vis, vengono scritti i dati nel file json e creata la pagina html con Leaflet
if product == "gprox":
    datiDaScriverenelJson = [NW_Latitude, NW_Longitude, SE_Latitude, SE_Longitude, timeIntervalStart, timeIntervalEnd, resolution, maxCloudCoverage, metersRadius, product]
    jsonBuilder(datiDaScriverenelJson, percentage_matrix, resultFolderPath)

if product == "stype":
    datiDaScriverenelJson = [NW_Latitude, NW_Longitude, SE_Latitude, SE_Longitude, timeIntervalStart, timeIntervalEnd, resolution, maxCloudCoverage, product]
    jsonBuilder(datiDaScriverenelJson, percentage_matrix, resultFolderPath)
    
if htmlFlag == "true":
    htmlLeafLetBuilder(NW_Longitude, NW_Latitude, SE_Longitude, SE_Latitude, resultFolderPath)


