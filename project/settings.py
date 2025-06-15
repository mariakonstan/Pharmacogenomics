# Project setup
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

FUZZY_IDENTITY_THRESHOLD = 0.9
EMAIL = "gene.scanner@gmail.com"

FASTA_DIR = BASE_DIR / "fasta"
FASTA_DIR.mkdir(parents=True, exist_ok=True)
FASTA_REPO = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome"
FASTA_NAME = "GRCh38_full_analysis_set_plus_decoy_hla"
FASTA_URLS = [
    f"{FASTA_REPO}/{FASTA_NAME}.fa",
    f"{FASTA_REPO}/{FASTA_NAME}.fa.fai",
    f"{FASTA_REPO}/{FASTA_NAME}.dict",
]

# Setup Input  
FASTA_FILE = f"{FASTA_DIR}/{FASTA_NAME}.fa" # Reference FASTA path
DATA_DIR = BASE_DIR / "data" # Input directory
GENES_JSON = DATA_DIR / "sequences/index.json" # Path to gene/allele JSON definition
SAMPLE_INDEX = DATA_DIR / "samples.index" # Sample index file (TSV)
POPULATIONS_JSON = DATA_DIR / "populations.json" # Path to population groups json

# Setup Output  
OUTPUT_DIR = BASE_DIR / "results" # Output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR_SAMPLES = OUTPUT_DIR / "samples" # Per sample output (json files)
OUTPUT_DIR_SAMPLES.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR_FIGURES = OUTPUT_DIR / "figures" # Plots (png files) 
OUTPUT_DIR_FIGURES.mkdir(parents=True, exist_ok=True)

# Setup logging  
LOG_DIR = OUTPUT_DIR / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOG_DIR / "scan.log"

# Run settings
RUN_FROM = 0 # Continue from iteration
RUN_LIMIT = 0 # Limit number of samples to scan (for debug)
RUN_PARALLEL = False # Run processing in parallel
RUN_PARALLEL_WORKERS = 4 # Number of parallel threads (if RUN_PARALLEL is True)
