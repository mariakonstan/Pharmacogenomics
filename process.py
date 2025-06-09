import project.settings as settings
from project.core import SampleScanner


def main():
    scanner = SampleScanner()
    scanner.genes_json = settings.GENES_JSON
    scanner.index_file = settings.SAMPLE_INDEX
    scanner.output_dir = settings.OUTPUT_DIR
    scanner.output_dir_samples = settings.OUTPUT_DIR_SAMPLES
    scanner.output_dir_figures = settings.OUTPUT_DIR_FIGURES
    scanner.load_gene_data()
    scanner.export_summaries()
    # scanner.plot_figures()
    
    
if __name__ == "__main__":
    main()