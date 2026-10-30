"""Download and unpack GEO GSE221553 (CD45-sorted scRNA-seq, pre/post TIL-ACT melanoma)."""
import tarfile
import urllib.request
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221553/suppl/GSE221553_RAW.tar"
RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
TAR = RAW / "GSE221553_RAW.tar"
OUT = RAW / "GSE221553"

RAW.mkdir(parents=True, exist_ok=True)
if not TAR.exists():
    print(f"Downloading {URL} (~165 MB)")
    urllib.request.urlretrieve(URL, TAR)
OUT.mkdir(exist_ok=True)
with tarfile.open(TAR) as t:
    t.extractall(OUT)
print(f"Unpacked {len(list(OUT.glob('*')))} files to {OUT}")
