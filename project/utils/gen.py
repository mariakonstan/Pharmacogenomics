import time
import pysam
import logging
import subprocess

from rapidfuzz import fuzz
from Bio import Entrez, SeqIO

from project.settings import EMAIL, FASTA_DIR, FASTA_URLS


Entrez.email = EMAIL
logger = logging.getLogger(__name__)


def hamming_identity(s1, s2):
    """
    Returns percent identity between two equal-length strings.
    """
    if len(s1) != len(s2):
        return 0.0
    matches = sum(c1 == c2 for c1, c2 in zip(s1, s2))
    return matches / len(s1)


def fuzz_similarity(s1, s2):
    """
    Returns rapidfuzz ratio between two non equal-length strings.
    """
    return fuzz.ratio(s1, s2) / 100
    

def clean_sequence(seq: str):
    """
    Removes line-breaks from sequence.
    """
    return seq.replace("\n", "")


def get_gene_sequence(accession, start=0, end=0, strand="", db="nucleotide"):
    """
    Returns Gene sequence from NCBI.
    """
    with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
        ix_start = start - 1 # Convert to 0-based indexing for Python slicing
        ix_end = end  # Python slices are exclusive at end
        seq = record.seq[ix_start:ix_end]
        return str(seq).upper()


def load_gene_sequence(file_path):
    """
    Returns Gene sequence from local data file.
    """
    with open(file_path, 'r') as f:
        return clean_sequence(f.read()).upper()

    
def download_fasta():
    """
    Downloads necessary input files.
    """
    for url in FASTA_URLS:
        filename = url.split("/")[-1]
        destination = FASTA_DIR / filename
        if not destination.exists():
            logger.info(f"Downloading {url} to {destination}")
            subprocess.run(["wget", "-O", str(destination), url], check=True)
        else:
            logger.info(f"File already exists: {destination}")    
   

def convert_ftp_to_https(url: str) -> str:
    if url.startswith("ftp://"):
        return url.replace("ftp://", "https://", 1)
    return url
   
   
def safe_open_alignment(url, ref_fasta, mode="rc", retries=3, delay=2):
    try:
        return pysam.AlignmentFile(url, mode=mode, reference_filename=ref_fasta)
    except Exception as e:
        if retries > 0:
            logger.warning(f"{url} > connection busy: Retrying... (remaining tries: {retries})")
            time.sleep(delay)
            return safe_open_alignment(url, ref_fasta, mode, retries-1, delay)
        