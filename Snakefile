"""
Author: Ivan Duque Aldana
iaduquea@unal.edu.co
Run: snakemake  -s Snakefile
"""

SAMPLES, = glob_wildcards("species1/{sample}.fasta")
LEN_SAMPLES = len(SAMPLES)

SAMPLES2, = glob_wildcards("species2/{sample}.fasta")

rule final:
  input: "clustering/0.7-identitygene.clstr",
            "clustering/0.7-identitygene",  
            "matrix/0.7_BIN_gene.npy",
            "matrix/0.7_FREQ_gene.npy",
            "matrix/0.7_PERC_gene.npy",
            "clustering/0.7-identityinter.clstr",
            "clustering/0.7-identityinter",
            "matrix/0.7_BIN_inter.npy",
            "matrix/0.7_FREQ_inter.npy",
            "matrix/0.7_PERC_inter.npy"     

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
    message: "-------Creating genome files required by bedtools----------"
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

rule concatenate_gene_sequences:
    input:  expand("Annotated/{sample}.genes", sample=SAMPLES),
        expand("Annotated/{sample}.genes", sample=SAMPLES2)
    output: "clustering/allGenes.fasta"
    message: "---------Concatenating genes sequences to clustering------------"
    shell: """cat {input} > {output}"""

rule concatenate_intergenic_sequences:
    input:  expand("intergenic/{sample}intergenic.fasta", sample=SAMPLES),
        expand("intergenic/{sample}intergenic.fasta", sample=SAMPLES2)
    output: "clustering/allintergenic.fasta"
    message: "---------Concatenating intergenic sequences to clustering------------"
    shell: """cat {input} > {output}"""

rule clustering_cd_hit_iterative_Gene:  
    input: "clustering/allGenes.fasta"
    output: "clustering/0.95-identitygene.clstr",
            "clustering/0.95-identitygene",
            "dics/0.95repsgene.json",
            "dics/0.95seqsgene.json",  
            "matrix/0.95_BIN_gene.npy",
            "matrix/0.95_FREQ_gene.npy",
            "matrix/0.95_PERC_gene.npy",
    params: region = "gene" , samples_file = 'dics/samples.json'
    message: "Iterative Clustering using cd-hit-est for genes"
    shell: "python scripts/clustering_cdhit.py {input[0]} {params.region} {params.samples_file}"


rule clustering_cd_hit_iterative_Intergenic:  
    input: "clustering/allintergenic.fasta"
    output: "clustering/0.95-identityinter.clstr",
            "clustering/0.95-identityinter",
            "dics/0.95repsinter.json",
            "dics/0.95seqsinter.json",
            "matrix/0.95_BIN_inter.npy",
            "matrix/0.95_FREQ_inter.npy",
            "matrix/0.95_PERC_inter.npy"
    params: region = "inter" , samples_file = 'dics/samples.json'
    message: "Iterative Clustering using cd-hit-est for intergenic"
    shell: "python scripts/clustering_cdhit.py {input[0]} {params.region} {params.samples_file}"

rule clustering_mcl__Genic:  
    input: "clustering/0.95-identitygene",
           "dics/0.95repsgene.json",
           "dics/0.95seqsgene.json"
    output: "clustering/0.7-identitygene.clstr",
            "clustering/0.7-identitygene",
            "matrix/0.7_BIN_gene.npy",
            "matrix/0.7_FREQ_gene.npy",
            "matrix/0.7_PERC_gene.npy",
    params: region = "gene" , samples_file = 'dics/samples.json'
    message: "Clustering applying all vs all using Blast and mcl for genes"
    shell: "python scripts/clustering_mcl.py {input[1]} {input[2]} {input[0]} {params.region} {params.samples_file}"


rule clustering_mcl_Intergenic:  
    input: "clustering/0.95-identityinter",
           "dics/0.95repsinter.json",
           "dics/0.95seqsinter.json"
    output: "clustering/0.7-identityinter.clstr",
            "clustering/0.7-identityinter",
            "matrix/0.7_BIN_inter.npy",
            "matrix/0.7_FREQ_inter.npy",
            "matrix/0.7_PERC_inter.npy"
    params: region = "inter" , samples_file = 'dics/samples.json'
    message: "Clustering applying all vs all using Blast and mcl for intergenic"
    shell: "python scripts/clustering_mcl.py {input[1]} {input[2]} {input[0]} {params.region} {params.samples_file}"""


    
