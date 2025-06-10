from project.core import SampleScanner
from project.utils import download_fasta
import project.settings as settings 


def main():
    download_fasta()
    
    scanner = SampleScanner()
    scanner.genes_json = settings.GENES_JSON
    scanner.index_file = settings.SAMPLE_INDEX
    scanner.ref_fasta_file = settings.FASTA_FILE 
    scanner.output_dir = settings.OUTPUT_DIR
    scanner.output_dir_samples = settings.OUTPUT_DIR_SAMPLES
    scanner.start = settings.RUN_FROM
    scanner.limit = settings.RUN_LIMIT
    scanner.load_gene_data()
    
    if settings.RUN_PARALLEL:
        scanner.run_parallel(max_workers=settings.RUN_PARALLEL_WORKERS)
    else:
        scanner.run()

    
if __name__ == "__main__":
    main()