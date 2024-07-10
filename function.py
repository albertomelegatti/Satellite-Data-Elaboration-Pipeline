import numpy as np
import cv2
from PIL import Image
import datetime
import json
import os
import tarfile
import datetime

# Funzione per ottenere una matrice di valori binari a partire da un'immagine
def getValuesMatrix(inputImagePath):
    image = cv2.imread(inputImagePath)
    height, width, _ = image.shape

    # La matrice viene inizializzata con tutti zeri
    value_matrix = np.zeros((height, width), dtype=np.int32)

    for y in range(height):
        for x in range(width):
            pixelColor = image[y, x]
            # Verifica se il colore del pixel è nero (RGB: 0, 0, 0)
            if np.array_equal(pixelColor, [0, 0, 0]):
                value_matrix[y][x] = 1

    return value_matrix

# Funzione per ottenere una matrice di percentuali di pixel verdi in un raggio intorno a ciascun pixel 
def getPercentageMatrix(matrix, radius):
    height, width = matrix.shape
    percentageMatrix = np.zeros((height, width), dtype=np.float64)
    total_pixel = height * width
    lastProgress = -1

    for y in range(height):
        center_y = y
        for x in range(width):

                greenCounter = 0
                totalCells = 0
                center_x = x

                for i in range(max(0, y - radius), min(height, y + radius + 1)):
                    for j in range(max(0, x - radius), min(width, x + radius + 1)):
                        distance = np.sqrt((j - center_x)**2 + (i - center_y)**2)
                        if distance <= radius:
                            totalCells += 1
                            if matrix[i][j] == 1:
                                greenCounter += 1

                if totalCells == 0:
                    greenPercentage = 0
                else:
                    greenPercentage = (greenCounter / totalCells) * 100
                    percentageMatrix[y][x] = greenPercentage


        # Stampa della percentuale di completamento
        progress = (y * width + x) / total_pixel * 100
        progress_int = int(progress)
        if progress_int > lastProgress:
            print('Progress: {:.0f}%'.format(progress_int), end='\r')
        lastProgress = progress_int
            
    return percentageMatrix

# A partire dalla matrice dei valori percentuali viene creata un'immagine a 5 bande di colore
def generateImage(matrix, outputPath):
    # Creazione di un'immagine vuota con le dimensioni della matrice
    height, width = matrix.shape
    image = np.ones((height, width, 3), dtype=np.uint8) * 255

    # Definizione delle soglie per le categorie e dei colori corrispondenti
    tresholds = [20, 40, 60, 80, 100]
    colours = [(255, 0, 0), (255, 128, 0), (255, 255, 0), (128, 255, 0), (0, 255, 0)]

    lastProgress = -1

    # Sovrascrittura dei valori della matrice sull'immagine con i colori corrispondenti alle categorie
    totalPixel = height * width
    for y in range(height):
        for x in range(width):
            # Trova la categoria corrispondente per il valore della matrice
            value = matrix[y][x]
            category = 0
            for treshold in tresholds:
                if value <= treshold:
                    colour = colours[category]
                    break
                category += 1

            # Imposta il colore del pixel nell'immagine
            image[y, x] = colour

            # Stampa della percentuale di completamento
            progress = (y * width + x) / totalPixel * 100
            progress_int = int(progress)
            if progress_int > lastProgress:
                print('Progress: {:.0f}%'.format(progress_int), end='\r')
            lastProgress = progress_int

    # Salvataggio dell'immagine
    image = Image.fromarray(image)
    image.save(outputPath, format='PNG')

# Funzione per confrontare due date in formato stringa
def compareDates(date_str1, date_str2):
    # Conversione delle date da stringa a oggetti datetime
    date1 = datetime.datetime.strptime(date_str1, '%Y-%m-%d').date()
    date2 = datetime.datetime.strptime(date_str2, '%Y-%m-%d').date()
    
    # 1 :  date1 viene prima di date2
    # 2 :  date1 viene dopo date2
    # 0 :  date1 e date2 sono uguali
    if date1 < date2:
        return 1
    elif date1 > date2:
        return 2
    else:
        return 0

# Funzione per confrontare la correttezza delle coordinate   
def compareCoordinates(lat1, lon1, lat2, lon2):
    # Controllo valori coordinate
    if lat1 < -90 or lat1 > 90 or lon1 < -180 or lon1 > 180 or lat2 < -90 or lat2 > 90 or lon2 < -180 or lon2 > 180:
        return 1
    # Controllo dell'ordine delle coordinate
    if lat1 < lat2 or lon1 > lon2:
        return 2
    return 0

# Funzione per estrarre le immagini dal file tar contenente le immagini
def extractImagesFromTar(outputFolderPath):
    # Ottieni il percorso corrente

    # Percorso della cartella "results"
    resultsDir = outputFolderPath

    # Entra nella cartella "results"
    os.chdir(resultsDir)

    # Ottieni la lista delle cartelle all'interno di "results"
    subfolders = [folder for folder in os.listdir() if os.path.isdir(folder)]

    # Controlla se c'è esattamente una cartella all'interno di "results"
    if len(subfolders) != 1:
        print("La cartella 'results' deve contenere esattamente una sotto-cartella.")
        return

    # Percorso della sola sotto-cartella
    subfolder_path = os.path.join(resultsDir, subfolders[0])

    # Ottieni la lista dei file nella sotto-cartella
    filesInSubfolder = os.listdir(subfolder_path)

    # Controlla se c'è esattamente un file tar nella sotto-cartella
    tar_files = [file for file in filesInSubfolder if file.endswith(".tar")]
    if len(tar_files) != 1:
        print("La sotto-cartella deve contenere esattamente un file tar.")
        return

    # Percorso completo del file tar
    tarFilePath = os.path.join(subfolder_path, tar_files[0])

    # Creazione della cartella dove verranno estratte le immagini
    extractionFolder = os.path.join(resultsDir, "extracted_contents")
    os.makedirs(extractionFolder, exist_ok=True)

    # Estrai il contenuto del file tar nella cartella di estrazione
    with tarfile.open(tarFilePath, 'r') as tar:
        tar.extractall(extractionFolder)

# Funzione che riceve in input un colore in formato RGB e la scelta effettuata e restituisce un valore diverso a seconda del colore e della scelta effettuata
def rgbToString(rgb):
    
    colore = '{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])
    if colore == "ff0000":
        return "Bare_Soil" #rosso
    elif colore == "00ff00":
        return "Vegetation" #verde
    elif colore == "ffffff":
        return "Building" #bianco
    elif colore == "0080ff":
        return "Water" #blu

def matrix2Json(dataToWrite, percentageMatrix, inputPath):
    outputPath = os.path.join(inputPath, "result.json")
    height, width = percentageMatrix.shape
    totalPixel = height * width

    data = {}
    data.update({"Product": dataToWrite[9]})
    data.update({"NW Latitude": dataToWrite[0]})
    data.update({"NW Longitude": dataToWrite[1]})
    data.update({"SE Latitude": dataToWrite[2]})
    data.update({"SE Longitude": dataToWrite[3]})
    data.update({"Time Interval Start": dataToWrite[4]})
    data.update({"Time Interval End": dataToWrite[5]})
    data.update({"Resolution (meter / pixel)": dataToWrite[6]})
    data.update({"Max Cloud Coverage (%)": dataToWrite[7]})
    data.update({"Radius (meter)": dataToWrite[8]})
    data.update({"Creation Date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})

    lastProgress = -1

    for y in range(height):
        key = "row" + str(y + 1)
        row_values = []
        for x in range(width):
            row_values.append(str(round(percentageMatrix[y][x], 2)))

            # Stampa della percentuale di completamento
            progress = (y * width + x) / totalPixel * 100
            progress_int = int(progress)
            if progress_int > lastProgress:
                print('Progress: {:.0f}%'.format(progress_int), end='\r')
            lastProgress = progress_int

        data[key] = row_values

    with open(outputPath, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# Funzione che scrive un file JSON a partire da un'immagine scaricata
def jsonBuilder(dataToWrite, percentageMtrix, inputPath):

    print("Writing Json file...")

    outputPath = os.path.join(inputPath, "result.json")

    if dataToWrite[len(dataToWrite) -1] == "gprox":
        matrix2Json(dataToWrite, percentageMtrix, inputPath)
        return 0

    elif dataToWrite[len(dataToWrite) -1] == "stemp":
        inputPath = os.path.join(inputPath, "extracted_contents","default.tif")
    
    else: 
        inputPath = os.path.join(inputPath, "extracted_contents", "convertedDefault.tif")

    im = cv2.imread(inputPath, cv2.IMREAD_UNCHANGED)
    im_data_rgb = cv2.cvtColor(im, cv2.COLOR_BGR2RGB)

    height, width = im_data_rgb.shape[:2]
    totalPixel = height * width

    # Viene creato un dizionario contenente le informazioni che andranno scritte nel file JSON
    dati = {"Product": dataToWrite[8]}
    dati.update({"NW Latitude": dataToWrite[0]})
    dati.update({"NW Longitude": dataToWrite[1]})
    dati.update({"SE Latitude": dataToWrite[2]})
    dati.update({"SE Longitude": dataToWrite[3]})
    dati.update({"Time Interval Start": dataToWrite[4]})
    dati.update({"Time Interval End": dataToWrite[5]})
    dati.update({"Resolution (meter / pixel)": dataToWrite[6]})
    dati.update({"Max Cloud Coverage (%)": dataToWrite[7]})
    dati.update({"Creation Date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})
    data = {}
    data.update(dati)

    lastProgress = -1

    # Lettura dell'immagine pixel per pixel e conversione in esadecimale
    for i in range(height):
        key = "row" + str(i + 1)
        row_values = []
        for j in range(width):
            pixel = im_data_rgb[i, j]
            row_values.append(rgbToString(pixel))

            # Stampa della percentuale di completamento
            progress = (i * width + j) / totalPixel * 100
            progress_int = int(progress)
            if progress_int > lastProgress:
                print('Progresso: {:.0f}%'.format(progress_int), end='\r')
            lastProgress = progress_int
            
        data[key] = row_values

    with open(outputPath, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# Funzione per convertire un'immagine in 4 colori ed eliminare i colori intermedi
def convertImageTo4Colors(inputPath, outputPath):

    print("Image Conversion...")

    # Definire i colori di destinazione
    colors = {
        "red": (255, 0, 0),
        "green": (0, 255, 0),
        "blue": (0, 128, 255),
        "white": (255, 255, 255)
    }

   # Funzione per trovare il colore più vicino
    def closestColor(pixel, colors):
        min_dist = float('inf')
        closest_color = None
        for color_value in colors.values():
            dist = np.linalg.norm(np.array(pixel) - np.array(color_value))
            if dist < min_dist:
                min_dist = dist
                closest_color = color_value
        return closest_color

    try:
        # Aprire l'immagine
        img = Image.open(inputPath)
        img = img.convert("RGB")  # Assicurarsi che l'immagine sia in RGB

        # Convertire l'immagine in un array numpy
        img_array = np.array(img)

        totalPixel = img_array.shape[0] * img_array.shape[1]
        lastProgress = -1

        # Approssimare ogni pixel al colore più vicino
        for i in range(img_array.shape[0]):
            for j in range(img_array.shape[1]):
                img_array[i, j] = closestColor(img_array[i, j], colors)

                
                # Stampa della percentuale di completamento
                progress = (i * img_array.shape[1] + j) / totalPixel * 100
                progress_int = int(progress)
                if progress_int > lastProgress:
                    print('Progress: {:.0f}%'.format(progress_int), end='\r')
                lastProgress = progress_int

        # Convertire l'array numpy di nuovo in un'immagine
        newImg = Image.fromarray(img_array.astype('uint8'), 'RGB')

        # Salvare l'immagine in formato TIFF
        newImg.save(outputPath, format='TIFF')

    except Exception as e:
        print(f"Conversion error: {e}")


# Funzione che converte un colore esadecimale in una temperatura in un range da 0 a 50
'''
def map_color_to_temperature(hex_color):
  
  # Rimuove il cancelletto iniziale dal colore esadecimale
  if hex_color.startswith('#'):
    hex_color = hex_color[1:]

  # Converte il colore esadecimale in un array NumPy di valori RGB
  rgb_color = np.array([int(hex_color[i:i+2], 16) for i in range(0, 6, 2)])

  # Calcola la temperatura
  temperature = (rgb_color[0] + rgb_color[1] + rgb_color[2]) / 3
  temperature = (temperature - 255) / -5.1
  temperature = np.clip(temperature, 0, 50)

  return temperature
'''

# Funzione per calcolare il punto medio tra due coordinate
def midPoint(longitudeNW, latitudeNW, longitudeSE, latitudeSE):
    return (latitudeNW + latitudeSE) / 2, (longitudeNW + longitudeSE) / 2

# Funzione per convertire un'immagine TIFF in formato JPG
def tiff2Jpg(inputPath, outputPath):
    try:
        # Apri l'immagine TIFF
        with Image.open(inputPath) as img:
            # Converti l'immagine in RGB (necessario per salvare come JPG)
            img = img.convert("RGB")
            # Salva l'immagine in formato JPG
            img.save(outputPath, "JPEG")
    except Exception as e:
        print(f"Conversion error: {e}")

# Funzione per creare un file HTML per visualizzare un'immagine su Leaflet
def htmlLeafLetBuilder(longitudeNW, latitudeNW, longitudeSE, latitudeSE, outputFolder):
    print("Creating HTML file...")
    midpoint_lat, midpoint_lon = midPoint(longitudeNW, latitudeNW, longitudeSE, latitudeSE)
    currentDir = outputFolder
    output_path = os.path.join(currentDir, "map.html")
    imagePath = os.path.join(currentDir, "extracted_contents", "default.jpg").replace("\\", "/")

    htmlCode = f'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Visualizzazione immagine su Leaflet</title>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <!-- Includi Leaflet CSS -->
        <link rel="stylesheet" href="https://unpkg.com/leaflet/dist/leaflet.css" />
        <!-- Stili per la mappa -->
        <style>
            html, body {{
                height: 100%;
                margin: 0;
                padding: 0;
            }}
            #map {{
                width: 100%;
                height: 100%;
            }}
            #slider {{
                position: absolute;
                top: 10px;
                left: 10px;
                z-index: 1000;
                background: white;
                padding: 10px;
                border-radius: 5px;
                box-shadow: 0 0 15px rgba(0,0,0,0.2);
            }}
        </style>
    </head>
    <body>
    
    <!-- Elemento div per la mappa -->
    <div id="map"></div>
    <!-- Slider per la trasparenza -->
    <div id="slider">
        <label for="opacityRange">Trasparenza: </label>
        <input type="range" id="opacityRange" min="0" max="1" step="0.01" value="0.4">
    </div>
    
    <!-- Includi Leaflet JS -->
    <script src="https://unpkg.com/leaflet/dist/leaflet.js"></script>
    
    <script>
        // Coordinate dei vertici della bounding box
        var bounds = [
            [{latitudeNW}, {longitudeNW}], // Nord-ovest
            [{latitudeSE}, {longitudeSE}]  // Sud-est
        ];
    
        // URL dell'immagine che vuoi visualizzare
        var imageUrl = '{imagePath}';
    
        // Crea la mappa Leaflet centrata sul punto medio
        var map = L.map('map').setView([{midpoint_lat}, {midpoint_lon}], 13);
    
        // Aggiungi l'overlay dell'immagine alla mappa
        var imageOverlay = L.imageOverlay(imageUrl, bounds, {{
            opacity: 0.4
        }}).addTo(map);
    
        // Aggiungi un layer di base (opzionale)
        L.tileLayer('https://{{s}}.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png', {{
            maxZoom: 19,
            attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }}).addTo(map);
    
        // Zoom sulla bounding box dell'immagine
        map.fitBounds(bounds);
    
        // Funzione per aggiornare la trasparenza dell'overlay dell'immagine
        document.getElementById('opacityRange').addEventListener('input', function(e) {{
            var opacity = e.target.value;
            imageOverlay.setOpacity(opacity);
        }});
    </script>
    
    </body>
    </html>
    '''

    # Scrivi il codice HTML nel file
    with open(output_path, 'w') as html_file:
        html_file.write(htmlCode)