import os
import requests
import zipfile
import shutil

URL = "https://zenodo.org/records/18302255/files/maps.zip"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ZIP_PATH = os.path.join(SCRIPT_DIR, "maps.zip")
EXTRACT_DIR = SCRIPT_DIR

def download_and_extract():
    print(f"Downloading from {URL}...")
    try:
        response = requests.get(URL, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0

        with open(ZIP_PATH, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    percent = (downloaded / total_size) * 100
                    print(f"Progress: {percent:.1f}% ({downloaded/(1024*1024):.1f}/{total_size/(1024*1024):.1f} MB)", end='\r')
                
        print(f"\nDownload complete. File saved as {ZIP_PATH}")
        
        print(f"Extracting {ZIP_PATH} to {EXTRACT_DIR}...")
        with zipfile.ZipFile(ZIP_PATH, 'r') as zip_ref:
            
            # Extract
            zip_ref.extractall(EXTRACT_DIR)
            
        print("Extraction complete.")
        
        # Cleanup zip file
        os.remove(ZIP_PATH)
        print(f"Removed {ZIP_PATH}")
        
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    download_and_extract()
