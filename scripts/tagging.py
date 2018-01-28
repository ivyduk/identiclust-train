from tag_fasta import tag_fasta

a = tag_fasta(snakemake.input.fastaSp1, snakemake.params.sp1, snakemake.params.start1)
b = tag_fasta(snakemake.input.fastaSp2, snakemake.params.sp2, snakemake.params.start2)
