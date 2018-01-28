# This should be placed in the Snakefile.

##-----------------------------------------------##
## Working directory                             ##
## Adapt to your needs                           ##
##-----------------------------------------------##

BASE_DIR = os.getcwd()
print (BASE_DIR)
##--------------------------------------------------------------------------------------##
## The list of samples to be processed
##--------------------------------------------------------------------------------------##
SAMPLES, = glob_wildcards("species1/{smp}.fasta")
LEN_SAMPLES = len(SAMPLES)

SAMPLES2, = glob_wildcards("species2/{smp}.fasta")

rule tag_samples:
    input:
        fastaSp1=expand("species1/{sample}.fasta", sample=SAMPLES),
        fastaSp2=expand("species2/{sample}.fasta", sample=SAMPLES2),
    output:
        'dics/SP1.json',
        'dics/SP2.json'
    message: "----Tagging samples of species folders-----"
    log:
        "logs/LOG.log"
    params: sp1 = "SP1", start1 = 0, sp2 = "SP2", start2 = LEN_SAMPLES
    script:
        "scripts/tagging.py"
