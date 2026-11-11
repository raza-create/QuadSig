import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_RAW       = os.path.join(PROJECT_ROOT, "data", "raw")
DATA_FASTA     = os.path.join(PROJECT_ROOT, "data", "fasta")
DATA_FASTA_MUT = os.path.join(PROJECT_ROOT, "data", "fasta", "mutants")
DATA_STRUCT    = os.path.join(PROJECT_ROOT, "data", "structures")
DATA_STRUCT_WT = os.path.join(PROJECT_ROOT, "data", "structures", "wildtype")
DATA_STRUCT_MUT= os.path.join(PROJECT_ROOT, "data", "structures", "mutants")
DATA_INTER     = os.path.join(PROJECT_ROOT, "data", "intermediate")
RESULTS_DIR    = os.path.join(PROJECT_ROOT, "results")
FIGURES_DIR    = os.path.join(PROJECT_ROOT, "figures")
TABLES_DIR     = os.path.join(PROJECT_ROOT, "tables")
LOGS_DIR       = os.path.join(PROJECT_ROOT, "logs")

UNIPROT_ID     = "P04637"
GENE_NAME      = "TP53"
SEQUENCE_LEN   = 393

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
NUM_MUTANTS = SEQUENCE_LEN * (len(AMINO_ACIDS) - 1)

DOMAINS = {
    "TAD": (1,   67,  "N-terminal Transactivation Domain"),
    "PRR": (68,  98,  "Proline-Rich Region"),
    "DBD": (102, 292, "DNA-Binding Domain"),
    "TET": (323, 356, "Tetramerisation Domain"),
    "CTD": (364, 393, "C-terminal Regulatory Domain"),
}

ZINC_RESIDUES = {176: "C", 179: "H", 238: "C", 242: "C"}

GOF_MUTATIONS = {
    "R175H": "Gains interaction with RAS effectors",
    "R248W": "DNA contact residue; GOF via MRE11 interaction",
    "R248Q": "Same position as R248W; different GOF mechanism",
    "R249S": "Hepatocellular carcinoma hotspot; gains interaction with CDK4",
    "R273H": "DNA contact residue; stabilises inactive conformation",
    "R273C": "Same position as R273H; GOF via SRC kinase",
    "G245S": "Structural mutation; GOF via interaction with NF-Y",
}

AA_CHARGE = {"A":0,"C":0,"D":-1,"E":-1,"F":0,"G":0,"H":0.1,"I":0,"K":1,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":1,"S":0,"T":0,"V":0,"W":0,"Y":0}
AA_MW = {"A":89,"C":121,"D":133,"E":147,"F":165,"G":75,"H":155,"I":131,"K":146,"L":131,"M":149,"N":132,"P":115,"Q":146,"R":174,"S":105,"T":119,"V":117,"W":204,"Y":181}
AA_HYDRO = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}

UNIPROT_API  = "https://rest.uniprot.org/uniprotkb"
ENSEMBL_API  = "https://rest.ensembl.org"
GNOMAD_API   = "https://gnomad.broadinstitute.org/api"
CLINVAR_API  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PATHOGENIC_THRESHOLD = 0.66
BENIGN_THRESHOLD     = 0.34
