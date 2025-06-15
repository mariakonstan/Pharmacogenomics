from enum import Enum


class SampleData:
    def __init__(self, name, url, pop=""):
        self.name = name
        self.url = url
        self.pop = pop


class Match:
    def __init__(self, sample_id, sample_pop, gene_name, match_type, score=None):
        """
        Match holds sample matching type and score.
        Args:
            sample_id (str): sample name or identifier
            sample_population (str): sample population 
            gene_name (str): gene/allele name or identifier
            match_type (str): characterization of matching (e.g. gene, gene_fuzzy, allele)
            score (float): 1.0 for complete match (< 1.0 fuzzy match)
        """
        self.sample_id = sample_id
        self.sample_population = sample_pop
        self.gene_name = gene_name 
        self.match_type = match_type  
        self.score = score  # could be frequency %, match % (max: 1.0)

    def to_dict(self):
        return {
            "sample_id": self.sample_id,
            "sample_population": self.sample_population,
            "gene_name": self.gene_name,
            "match_type": self.match_type,
            "score": self.score
        }


class GeneType(Enum):
    NONE = ""
    GENE = "GENE"
    ALLELE = "ALLELE"


class Gene:
    def __init__(self, name, chrom="", start=0, end=0):
        """
        GeneData holds genomic region info.
        Args:
            name (str): gene name or identifier
            chrom (str): chromosome name (e.g. 'chr17')
            start (int): 1-based start position
            end (int): end position
        """
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.type = GeneType.NONE
        self.pos = None # if defined it's a gene allele
        self.seq = "" # nucleotide sequence (len should be end-start)
        self.tag = "" # some label to group gene-data
        self.attr = "" # characterization depending on drug response
        self.comments = "" 
        # Per sample evaluation 
        self.match = None # sample match
        
        
class GeneGroup:
    def __init__(self, tag, gene:Gene, alleles=[]):
        self.tag = tag
        self.gene = gene
        self.alleles = alleles
