import os
import requests
import zipfile
import shutil
import subprocess
import sys

URL = "https://zenodo.org/records/18302255/files/maps.zip"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ZIP_PATH = os.path.join(SCRIPT_DIR, "maps.zip")
EXTRACT_DIR = SCRIPT_DIR

def download_and_extract():
    print(f"Downloading from {URL}...")
    
    downloaded_via_system = False
    
    # Method 1: Try system curl (often more robust on Linux/Mac)
    if shutil.which("curl"):
        print("Detected 'curl'. Attempting download via system command...")
        try:
            subprocess.run(
                ["curl", "-L", "-o", ZIP_PATH, URL], 
                check=True
            )
            downloaded_via_system = True
            print("\nCurl download complete.")
        except subprocess.CalledProcessError:
            print("Curl download failed. Falling back to Python requests.")
        except Exception as e:
            print(f"Curl execution error: {e}. Falling back to Python requests.")

    # Method 2: Python requests (Windows fallback or if curl fails)
    if not downloaded_via_system:
        try:
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            
            response = requests.get(URL, headers=headers, stream=True, timeout=30)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0

            with open(ZIP_PATH, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"Progress: {percent:.1f}% ({downloaded/(1024*1024):.1f}/{total_size/(1024*1024):.1f} MB)", end='\r', flush=True)
            
            print(f"\nDownload complete. File saved as {ZIP_PATH}")
            
        except Exception as e:
            print(f"An error occurred during requests download: {e}")
            return

    print(f"Extracting {ZIP_PATH} to {EXTRACT_DIR}...")
    with zipfile.ZipFile(ZIP_PATH, 'r') as zip_ref:
        
        # Extract
        zip_ref.extractall(EXTRACT_DIR)
        
    print("Extraction complete.")
    
    os.remove(ZIP_PATH)
    print(f"Removed {ZIP_PATH}")

if __name__ == "__main__":
    download_and_extract()
