"""
Author: Ivan Duque Aldana
iaduquea@unal.edu.co
Run: snakemake   -s Snakefile
"""

##--------------------------------------------------------------------------------------##
## The list of samples to be processed
##--------------------------------------------------------------------------------------##
SAMPLES, = glob_wildcards("species1/{sample}.fasta")
LEN_SAMPLES = len(SAMPLES)

SAMPLES2, = glob_wildcards("species2/{sample}.fasta")

rule final:
  input: expand("Annotated/{sample}.gff", sample=SAMPLES),
         expand("Annotated/{sample}.genes", sample=SAMPLES)

rule tag_samples:
    input:
        fastaSp1=expand("species1/{sample}.fasta", sample=SAMPLES),
        fastaSp2=expand("species2/{sample}.fasta", sample=SAMPLES2),
    output:
        'dics/SP1.json',
        'dics/SP2.json'
    message: "---- Tagging samples of species folders -----"
    params: sp1 = "SP1", start1 = 0, sp2 = "SP2", start2 = LEN_SAMPLES
    script:
        "scripts/tagging.py"

rule gene_annotation:
    input: "species1/{sample}.fasta"
    output:
        gff = 'Annotated/{sample}.gff', genes = 'Annotated/{sample}.genes'
    message: "---- Predicting genes with prodigal -----"
    shell: """prodigal -i {input} -m -f gff -o {output.gff} -d {output.genes}"""
