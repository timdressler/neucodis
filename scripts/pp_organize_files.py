import os
import shutil
import re

# Setze den Pfad zu deinem Ordner mit den Dateien
source_folder = "DEIN_PFAD_ZUM_QUELLORDNER"
# Setze den Pfad zum Zielordner, in dem die neuen Ordner erstellt werden sollen
destination_folder = "DEIN_PFAD_ZUM_ZIELORDNER"

# Erstelle eine Liste aller Dateien im Quellordner
files = os.listdir(source_folder)

# Durchlaufe alle Dateien im Ordner
for file_name in files:
    # Überprüfe, ob der Dateiname das Muster "av_PXXX_" erfüllt
    match = re.match(r"av_(P\d{3})_", file_name)
    
    if match:
        # Extrahiere den Ordnernamen (z.B. "P120")
        folder_name = match.group(1)
        target_folder = os.path.join(destination_folder, folder_name)
        
        # Erstelle den Zielordner, falls er nicht existiert
        os.makedirs(target_folder, exist_ok=True)
        
        # Voller Pfad zur Datei
        source_file = os.path.join(source_folder, file_name)
        target_file = os.path.join(target_folder, file_name)
        
        # Verschiebe die Datei in den entsprechenden Ordner im Zielverzeichnis
        shutil.move(source_file, target_file)

print("Dateien wurden erfolgreich sortiert und verschoben.")
