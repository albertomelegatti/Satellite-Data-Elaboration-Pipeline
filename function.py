import numpy as np
import cv2
from PIL import Image
import datetime
import json
import os
import tarfile
import datetime

# Funzione per ottenere una matrice di valori binari a partire da un'immagine
def getValuesMatrix(percorso_immagine_input):
    immagine = cv2.imread(percorso_immagine_input)
    height, width, _ = immagine.shape

    # La matrice viene inizializzata con tutti zeri
    value_matrix = np.zeros((height, width), dtype=np.int32)

    for y in range(height):
        for x in range(width):
            colore_pixel = immagine[y, x]
            # Verifica se il colore del pixel è nero (RGB: 0, 0, 0)
            if np.array_equal(colore_pixel, [0, 0, 0]):
                value_matrix[y][x] = 1

    return value_matrix

# Funzione per ottenere una matrice di percentuali di pixel verdi in un raggio intorno a ciascun pixel 
def getPercentageMatrix(matrice, raggio):
    altezza, larghezza = matrice.shape
    matrice_percentuali = np.zeros((altezza, larghezza), dtype=np.float64)
    total_pixel = altezza * larghezza
    ultimo_progresso = -1

    for y in range(altezza):
        centro_y = y
        for x in range(larghezza):

                conta_verdi = 0
                totale_celle = 0
                centro_x = x

                for i in range(max(0, y - raggio), min(altezza, y + raggio + 1)):
                    for j in range(max(0, x - raggio), min(larghezza, x + raggio + 1)):
                        distanza = np.sqrt((j - centro_x)**2 + (i - centro_y)**2)
                        if distanza <= raggio:
                            totale_celle += 1
                            if matrice[i][j] == 1:
                                conta_verdi += 1

                if totale_celle == 0:
                    percentuale_verdi = 0
                else:
                    percentuale_verdi = (conta_verdi / totale_celle) * 100
                    matrice_percentuali[y][x] = percentuale_verdi


        # Stampa della percentuale di completamento
        progresso = (y * larghezza + x) / total_pixel * 100
        progresso_int = int(progresso)
        if progresso_int > ultimo_progresso:
            print('Progresso: {:.0f}%'.format(progresso_int), end='\r')
        ultimo_progresso = progresso_int
            
    return matrice_percentuali

# A partire dalla matrice dei valori percentuali viene creata un'immagine a 5 bande di colore
def generateImage(matrice, percorso_output):
    # Creazione di un'immagine vuota con le dimensioni della matrice
    altezza, larghezza = matrice.shape
    immagine = np.ones((altezza, larghezza, 3), dtype=np.uint8) * 255

    # Definizione delle soglie per le categorie e dei colori corrispondenti
    soglie = [20, 40, 60, 80, 100]
    colori = [(255, 0, 0), (255, 128, 0), (255, 255, 0), (128, 255, 0), (0, 255, 0)]

    ultimo_progresso = -1

    # Sovrascrittura dei valori della matrice sull'immagine con i colori corrispondenti alle categorie
    totali_pixel = altezza * larghezza
    for y in range(altezza):
        for x in range(larghezza):
            # Trova la categoria corrispondente per il valore della matrice
            valore = matrice[y][x]
            categoria = 0
            for soglia in soglie:
                if valore <= soglia:
                    colore = colori[categoria]
                    break
                categoria += 1

            # Imposta il colore del pixel nell'immagine
            immagine[y, x] = colore

            # Stampa della percentuale di completamento
            progresso = (y * larghezza + x) / totali_pixel * 100
            progresso_int = int(progresso)
            if progresso_int > ultimo_progresso:
                print('Progresso: {:.0f}%'.format(progresso_int), end='\r')
            ultimo_progresso = progresso_int

    # Salvataggio dell'immagine
    immagine = Image.fromarray(immagine)
    immagine.save(percorso_output, format='PNG')

# Funzione per confrontare due date in formato stringa
def compare_dates(date_str1, date_str2):
    # Convert the string dates to datetime.date objects
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
def compare_coordinates(lat1, lon1, lat2, lon2):
    # Check if the coordinates are within the valid range
    if lat1 < -90 or lat1 > 90 or lon1 < -180 or lon1 > 180 or lat2 < -90 or lat2 > 90 or lon2 < -180 or lon2 > 180:
        return 1
    # Check if the coordinates are in the correct order
    if lat1 < lat2 or lon1 > lon2:
        return 2
    return 0

# Funzione per estrarre le immagini dal file tar contenente le immagini
def extractImagesFromTar(posizioneCartellaRisultato):
    # Ottieni il percorso corrente

    # Percorso della cartella "results"
    results_dir = posizioneCartellaRisultato

    # Entra nella cartella "results"
    os.chdir(results_dir)

    # Ottieni la lista delle cartelle all'interno di "results"
    subfolders = [folder for folder in os.listdir() if os.path.isdir(folder)]

    # Controlla se c'è esattamente una cartella all'interno di "results"
    if len(subfolders) != 1:
        print("La cartella 'results' deve contenere esattamente una sotto-cartella.")
        return

    # Percorso della sola sotto-cartella
    subfolder_path = os.path.join(results_dir, subfolders[0])

    # Ottieni la lista dei file nella sotto-cartella
    files_in_subfolder = os.listdir(subfolder_path)

    # Controlla se c'è esattamente un file tar nella sotto-cartella
    tar_files = [file for file in files_in_subfolder if file.endswith(".tar")]
    if len(tar_files) != 1:
        print("La sotto-cartella deve contenere esattamente un file tar.")
        return

    # Percorso completo del file tar
    tar_file_path = os.path.join(subfolder_path, tar_files[0])

    # Creazione della cartella dove verranno estratte le immagini
    extraction_folder = os.path.join(results_dir, "extracted_contents")
    os.makedirs(extraction_folder, exist_ok=True)

    # Estrai il contenuto del file tar nella cartella di estrazione
    with tarfile.open(tar_file_path, 'r') as tar:
        tar.extractall(extraction_folder)

    #print(f"Il contenuto del file tar è stato estratto in: {extraction_folder}")

# Funzione che riceve in input un colore in formato RGB e la scelta effettuata e restituisce un valore diverso a seconda del colore e della scelta effettuata
def rgb_to_hex(rgb, scelta):
    
    colore = '{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])
    if scelta == "stype":
        if colore == "ff0000":
            return "Bare_Soil" #rosso
        elif colore == "00ff00":
            return "Vegetation" #verde
        elif colore == "ffffff":
            return "Building" #bianco
        elif colore == "0080ff":
            return "Water" #blu
        
    elif scelta == "stemp":
        temperatura = map_color_to_temperature(colore)
        temperaturaStringa = str(round(temperatura, 2))
        return temperaturaStringa

def matrix2Json(datiDaScrivere, matricePercentuali, inputPath):
    outputPath = os.path.join(inputPath, "result.json")
    altezza, larghezza = matricePercentuali.shape
    totali_pixel = altezza * larghezza

    data = {}
    data.update({"Product": datiDaScrivere[9]})
    data.update({"NW Latitude": datiDaScrivere[0]})
    data.update({"NW Longitude": datiDaScrivere[1]})
    data.update({"SE Latitude": datiDaScrivere[2]})
    data.update({"SE Longitude": datiDaScrivere[3]})
    data.update({"Time Interval Start": datiDaScrivere[4]})
    data.update({"Time Interval End": datiDaScrivere[5]})
    data.update({"Resolution (meter / pixel)": datiDaScrivere[6]})
    data.update({"Max Cloud Coverage (%)": datiDaScrivere[7]})
    data.update({"Radius (meter)": datiDaScrivere[8]})
    data.update({"Creation Date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})

    ultimo_progresso = -1

    for y in range(altezza):
        key = "row" + str(y + 1)
        row_values = []
        for x in range(larghezza):
            row_values.append(str(round(matricePercentuali[y][x], 2)))

            # Stampa della percentuale di completamento
            progresso = (y * larghezza + x) / totali_pixel * 100
            progresso_int = int(progresso)
            if progresso_int > ultimo_progresso:
                print('Progresso: {:.0f}%'.format(progresso_int), end='\r')
            ultimo_progresso = progresso_int

        data[key] = row_values

    with open(outputPath, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# Funzione che scrive un file JSON a partire da un'immagine scaricata
def jsonBuilder(datiDaScrivere, matricePercentuali, inputPath):

    print("Scrittura del file JSON in corso...")

    outputPath = os.path.join(inputPath, "result.json")

    if datiDaScrivere[len(datiDaScrivere) -1] == "gprox":
        matrix2Json(datiDaScrivere, matricePercentuali, inputPath)
        return 0

    elif datiDaScrivere[len(datiDaScrivere) -1] == "stemp":
        inputPath = os.path.join(inputPath, "extracted_contents","default.tif")
    
    else: 
        inputPath = os.path.join(inputPath, "extracted_contents", "convertedDefault.tif")

    im = cv2.imread(inputPath, cv2.IMREAD_UNCHANGED)
    im_data_rgb = cv2.cvtColor(im, cv2.COLOR_BGR2RGB)

    height, width = im_data_rgb.shape[:2]
    totali_pixel = height * width

    # Viene creato un dizionario contenente le informazioni che andranno scritte nel file JSON
    dati = {"Product": datiDaScrivere[8]}
    dati.update({"NW Latitude": datiDaScrivere[0]})
    dati.update({"NW Longitude": datiDaScrivere[1]})
    dati.update({"SE Latitude": datiDaScrivere[2]})
    dati.update({"SE Longitude": datiDaScrivere[3]})
    dati.update({"Time Interval Start": datiDaScrivere[4]})
    dati.update({"Time Interval End": datiDaScrivere[5]})
    dati.update({"Resolution (meter / pixel)": datiDaScrivere[6]})
    dati.update({"Max Cloud Coverage (%)": datiDaScrivere[7]})
    dati.update({"Creation Date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})
    data = {}
    data.update(dati)

    ultimo_progresso = -1

    # Lettura dell'immagine pixel per pixel e conversione in esadecimale
    for i in range(height):
        key = "row" + str(i + 1)
        row_values = []
        for j in range(width):
            pixel = im_data_rgb[i, j]
            row_values.append(rgb_to_hex(pixel, datiDaScrivere[8]))

            # Stampa della percentuale di completamento
            progresso = (i * width + j) / totali_pixel * 100
            progresso_int = int(progresso)
            if progresso_int > ultimo_progresso:
                print('Progresso: {:.0f}%'.format(progresso_int), end='\r')
            ultimo_progresso = progresso_int
            
        data[key] = row_values

    with open(outputPath, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# Funzione per convertire un'immagine in 4 colori
def convert_image_to_4_colors(inputPath, outputPath):

    print("Conversione dell 'immagine in corso...")

    # Definire i colori di destinazione
    colors = {
        "red": (255, 0, 0),
        "green": (0, 255, 0),
        "blue": (0, 128, 255),
        "white": (255, 255, 255)
    }

   # Funzione per trovare il colore più vicino
    def closest_color(pixel, colors):
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

        totali_pixel = img_array.shape[0] * img_array.shape[1]
        ultimo_progresso = -1

        # Approssimare ogni pixel al colore più vicino
        for i in range(img_array.shape[0]):
            for j in range(img_array.shape[1]):
                img_array[i, j] = closest_color(img_array[i, j], colors)

                
                # Stampa della percentuale di completamento
                progresso = (i * img_array.shape[1] + j) / totali_pixel * 100
                progresso_int = int(progresso)
                if progresso_int > ultimo_progresso:
                    print('Progresso: {:.0f}%'.format(progresso_int), end='\r')
                ultimo_progresso = progresso_int

        # Convertire l'array numpy di nuovo in un'immagine
        new_img = Image.fromarray(img_array.astype('uint8'), 'RGB')

        # Salvare l'immagine in formato TIFF
        new_img.save(outputPath, format='TIFF')
        #print(f"Immagine salvata correttamente come {outputPath}")
    except Exception as e:
        print(f"Errore durante la conversione: {e}")

# Funzione che converte un colore esadecimale in una temperatura in un range da 0 a 50
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

# Funzione per calcolare il punto medio tra due coordinate
def midpoint(longitudineNO, latitudineNO, longitudineSE, latitudineSE):
    return (latitudineNO + latitudineSE) / 2, (longitudineNO + longitudineSE) / 2

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
        print(f"Errore durante la conversione dell'immagine: {e}")

# Funzione per creare un file HTML per visualizzare un'immagine su Leaflet
def htmlLeafLetBuilder(longitudineNO, latitudineNO, longitudineSE, latitudineSE, outputFolder):
    midpoint_lat, midpoint_lon = midpoint(longitudineNO, latitudineNO, longitudineSE, latitudineSE)
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
            [{latitudineNO}, {longitudineNO}], // Nord-ovest
            [{latitudineSE}, {longitudineSE}]  // Sud-est
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

    with open(output_path, 'w') as html_file:
        html_file.write(htmlCode)