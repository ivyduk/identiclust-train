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
  input: expand("intergenic/{sample}intergenic.fasta", sample=SAMPLES),
        expand("intergenic/{sample}intergenic.fasta", sample=SAMPLES2)

rule tag_samples:
    input:
        fastaSp1=expand("species1/{sample}.fasta", sample=SAMPLES),
        fastaSp2=expand("species2/{sample}.fasta", sample=SAMPLES2),
    output:
        'dics/samples.json',
        expand("Genomes/{sample}.fasta", sample=SAMPLES),
        expand("Genomes/{sample}.fasta", sample=SAMPLES2)
    message: "---- Tagging samples of species folders -----"
    params: sp1 = "SP1", start1 = 0, sp2 = "SP2", start2 = LEN_SAMPLES
    script:
        "scripts/tagging.py"

rule gene_annotation:
    input: "Genomes/{sample}.fasta"
    output:
        gff = 'Annotated/{sample}.gff', genes = 'Annotated/{sample}.genes'
    message: "---- Predicting genes with prodigal -----"
    shell: """prodigal -i {input} -m -f gff -o {output.gff} -d {output.genes}"""

rule genome_file_creation:
    input: "Annotated/{sample}.gff"
    output:
        "Genomes/{sample}.genome"
    message: "---------Creating genome files------------"
    run: 
            import csv, re   
            f = open(input[0], "r")
            strings = re.findall(r'seqlen=(\d+);seqhdr="(\S+).*\"', f.read())
            with open(output[0], 'w') as csvfile:
                print ("writing genome file")
                writer = csv.writer(csvfile, lineterminator="\n", delimiter='\t')
                writer.writerow(['chrom','size'])
                for string in strings:
                    writer.writerow([string[1],string[0]])

rule intergenic_extraction_bed:
    input: "Genomes/{sample}.genome",
            "Annotated/{sample}.gff"
    output: "intergenic/{sample}.bed"
    message: "---------Creating bed files to obtain intergenic sequences------------"
    shell: """bedtools complement -i {input[1]} -g {input[0]} | awk '(($2 < $3) || (($2 == $3) && ($2 != 0)))' > {output}"""

rule intergenic_extraction_fasta:
    input:  "intergenic/{sample}.bed",
            "Genomes/{sample}.fasta"
    output:  "intergenic/{sample}intergenic.fasta",
    message: "---------Creating fasta files with intergenic sequences------------"
    shell: """bedtools getfasta -fi {input[1]} -bed {input[0]} -fo {output}"""

    
