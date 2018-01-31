from tag_fasta import tag_fasta
import json, shutil, os

a = tag_fasta(snakemake.input.fastaSp1, snakemake.params.sp1, snakemake.params.start1)
b = tag_fasta(snakemake.input.fastaSp2, snakemake.params.sp2, snakemake.params.start2)

b.update(a)

try:
    os.stat("Genomes")
except:
    os.mkdir("Genomes")

dest = 'Genomes/'
for f in snakemake.input.fastaSp1:
    shutil.move(f, dest)
for f in snakemake.input.fastaSp2:
    shutil.move(f, dest)

with open("dics/samples.json", 'w') as f:
        json_str = json.dumps(b, indent=4)
        f.write(json_str)
        f.close()
