"""
Author: Ivan Duque Aldana
iaduquea@unal.edu.co
Run: snakemake  -s Snakefile
"""

SAMPLES, = glob_wildcards("species1/{sample}.fasta")
LEN_SAMPLES = len(SAMPLES)

SAMPLES2, = glob_wildcards("species2/{sample}.fasta")
LEN_SAMPLES2 = len(SAMPLES2)

rule final:
  input: "models/finalized_model_gene.sav",
         "models/finalized_model_inter.sav",
         "models/finalized_model_both.sav"   

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
    shell: """prodigal -i {input} -m -f gff -o {output.gff} -d {output.genes} -c"""

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
    output: "clustering/0.8-identitygene.clstr",
            "clustering/0.8-identitygene",
            "dics/0.8clstrsgene.json",
            "dics/0.8seqsgene.json",  
            "matrix/0.8_BIN_gene.npy",
            "matrix/0.8_FREQ_gene.npy",
            "matrix/0.8_PERC_gene.npy",
    params: region = "gene" , samples_file = 'dics/samples.json'
    message: "Iterative Clustering using cd-hit-est for genes"
    shell: "python scripts/clustering_cdhit.py {input[0]} {params.region} {params.samples_file}"


rule clustering_cd_hit_iterative_Intergenic:  
    input: "clustering/allintergenic.fasta"
    output: "clustering/0.8-identityinter.clstr",
            "clustering/0.8-identityinter",
            "dics/0.8clstrsinter.json",
            "dics/0.8seqsinter.json",
            "matrix/0.8_BIN_inter.npy",
            "matrix/0.8_FREQ_inter.npy",
            "matrix/0.8_PERC_inter.npy"
    params: region = "inter" , samples_file = 'dics/samples.json'
    message: "Iterative Clustering using cd-hit-est for intergenic"
    shell: "python scripts/clustering_cdhit.py {input[0]} {params.region} {params.samples_file}"

rule clustering_mcl__Genic:  
    input: "clustering/0.8-identitygene",
           "dics/0.8clstrsgene.json",
           "dics/0.8seqsgene.json",
            "clustering/allGenes.fasta"           
    output: "clustering/0.7-identitygene.clstr",
            "clustering/0.7-identitygene",
            "matrix/0.7_BIN_gene.npy",
            "matrix/0.7_FREQ_gene.npy",
            "matrix/0.7_PERC_gene.npy",
    params: region = "gene" , samples_file = 'dics/samples.json'
    message: "Clustering applying all vs all using Blast and mcl for genes"
    shell: "python scripts/clustering_mcl.py {input[1]} {input[2]} {input[0]} {params.region} {params.samples_file} {input[3]}"


rule clustering_mcl_Intergenic:  
    input: "clustering/0.8-identityinter",
           "dics/0.8clstrsinter.json",
           "dics/0.8seqsinter.json",
           "clustering/allintergenic.fasta"
    output: "clustering/0.7-identityinter.clstr",
            "clustering/0.7-identityinter",
            "matrix/0.7_BIN_inter.npy",
            "matrix/0.7_FREQ_inter.npy",
            "matrix/0.7_PERC_inter.npy"
    params: region = "inter" , samples_file = 'dics/samples.json'
    message: "Clustering applying all vs all using Blast and mcl for intergenic"
    shell: "python scripts/clustering_mcl.py {input[1]} {input[2]} {input[0]} {params.region} {params.samples_file} {input[3]}"""

rule train_model_random_forest_gene:
    input: "matrix/0.7_PERC_gene.npy"
    output: "models/finalized_model_gene.sav"
    params: region = 'gene' ,samples_file = 'dics/samples.json', samples_SP2_number = LEN_SAMPLES2
    message: "Training Model for genic clusters using Random Forest"
    shell:"python scripts/train_model.py {input[0]} {params.region} {params.samples_file} {params.samples_SP2_number}"""


rule train_model_random_forest_inter:
    input: "matrix/0.7_PERC_inter.npy"
    output: "models/finalized_model_inter.sav"
    params: region = 'inter' ,samples_file = 'dics/samples.json', samples_SP2_number = LEN_SAMPLES2
    message: "Training Model for intergenic clusters using Random Forest"
    shell:"python scripts/train_model.py {input[0]} {params.region} {params.samples_file} {params.samples_SP2_number}"""

rule join_gene_and_intergenic_matrix:
    input: "matrix/0.7_PERC_gene.npy",
            "matrix/0.7_PERC_inter.npy"
    output: "matrix/0.7_PERC_both.npy"
    run: 
        import numpy as np
        TrainGene=np.load(input[0])
        TrainInter=np.load(input[1])
        TotalMatrixTrain = np.concatenate((TrainGene, TrainInter), axis=1)
        np.save("matrix/0.7_PERC_both", TotalMatrixTrain)

rule train_model_random_forest_both:
    input: "matrix/0.7_PERC_both.npy"
    output: "models/finalized_model_both.sav"
    params: region = 'both' ,samples_file = 'dics/samples.json', samples_SP2_number = LEN_SAMPLES2
    message: "Training Model for both types of clusters using Random Forest"
    shell:"python scripts/train_model.py {input[0]} {params.region} {params.samples_file} {params.samples_SP2_number}"""



    


    
