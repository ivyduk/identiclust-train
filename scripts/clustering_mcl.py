from clustering_cdhit import readJson, saveJson
import json
import os
import sys
import subprocess
import csv
import operator
from Bio import SeqIO
from clustering_cdhit import buildMatrix



def allvsallclustering(representatives, sequences, file_reps, region_type, samples_file): 
	
	samples_dic = readJson(samples_file)
	samplesnumber = len(samples_dic.keys())

	directoryallvsall = 'allVSall70'
	if not os.path.exists(directoryallvsall):
		os.makedirs(directoryallvsall)
	subprocess.call(['cp', file_reps, 'allVSall70/'])

	file_name = file_reps.split('/')[-1]
	print (file_name, file_reps)
	print ("creating blast database for Representatives comparison")
	subprocess.call(['makeblastdb', '-in', '%s/%s' % (directoryallvsall,file_name), '-parse_seqids', '-dbtype', 'nucl', '-out', '%s/%s%s' % (directoryallvsall,file_name,region_type) ])
	print ("Blast Database created for Sequences Representatives of previous clustering proccess")
	#execute blast for sequences against itself
	subprocess.call(['blastn', '-db',  '%s/%s%s' % (directoryallvsall,file_name,region_type), '-query', file_reps, '-perc_identity', '%s' % str(70),'-out', '%s/%sallblasted' % (directoryallvsall,region_type), '-outfmt', "6 qseqid sseqid bitscore qlen length pident qcovs" ])
	#dictionary with length of each sequence used in BLAST coomparison
	tuples_len = {}
	#dictionary that stores tubles as keys of two sequences and the value is the identity calculated
	tuple_dic = {}
	border="---8<------8<------8<------8<------8<---" + "\n"
	file = open("%s/mcl_input%s" % (directoryallvsall,region_type),"a")
	file.write(border)
	with open('%s/%sallblasted' % (directoryallvsall,region_type), 'r') as csvfile:
		blastedreader = csv.reader(csvfile, delimiter='\t')
		for row in blastedreader:
			tuples_len[row[0]]=row[3]
			tuples_len[row[1]]=row[4]
			file.write(row[0]+"\t"+row[1]+"\t"+row[5]+"\n")
			# example {(seq1, seq2) : '95.87'
			tuple_dic[(row[0],row[1])]=row[5]
			tuple_dic[(row[1],row[0])]=row[5]
	file.write(border)
	file.close()
	#Execute clustering with MCL
	print ("executing MCL algortihm")
	subprocess.call(['mcl', "%s/mcl_input%s" % (directoryallvsall,region_type), '-I', '1.4','--abc', '-o', "%s/mcl_out%s.all" % (directoryallvsall,region_type)])
	#Build clustering file output
	print ("Building clustering file output")
	file_mcl = open("%s/mcl_out%s.all" % (directoryallvsall,region_type),"rU")
	clusters_dic = {}
	for line in file_mcl.readlines():
		dic_temp = {}
		elements = []
		elements = line.split('\t')
		elements[-1] = elements[-1].strip()
		for i in elements:
			a=i
			dic_temp[a]=tuples_len[a]
		#extract representative by longest length 
		rep = max(dic_temp.items(), key=operator.itemgetter(1))[0]
		elements.remove(rep)
		for n,i in enumerate(elements):
			# if in blast results there is a identity measure among representative and sequence in cluster 
			if (rep,i) in tuple_dic.keys():
				elements[n]=(i,tuple_dic[(rep,i)])
			else:
				elements[n]=(i,'80.00')
		clusters_dic[rep]=elements
	with open("clustering/0.7-identity%s.clstr" %region_type, 'w') as f:
		json_str = json.dumps(clusters_dic, indent=4)
		f.write(json_str)
	print (' %s clusters identified with mcl method' % str(len(clusters_dic.keys())))
	
	
	#Updating files
	cluster = 0
	sequences = readJson(sequences)
	total_sequences = []
	new_representatives = {}
	new_representatives
	prev_representatives = readJson(representatives)
	
	#Update Representatives file
	for rep, results in clusters_dic.items():
		if rep not in total_sequences:
			total_sequences.append(rep)
		seqs = results
		representative = rep
		#sequences stored in previous cluster of representatvive
		prev_key_sequences = prev_representatives[representative]['sequences'] 
		key_sequences = prev_key_sequences 
		new_representatives[representative] = {}
		new_representatives[representative]["cluster_id"] = str(cluster)

		#update cluster info in representative
		sequences[rep]['cluster_id'] = str(cluster)
		sequences[rep]['identity_to_rep'] = '100.00'

		if not results:
			for r in prev_key_sequences:
				sequences[r]['cluster_id'] = str(cluster)

		for s in seqs:	
			if s[0] not in total_sequences:
				total_sequences.append(s[0])

			if s[0] not in key_sequences:
				key_sequences.append(s[0])
			#update sequences in previous representatives
			for prev_seq in prev_representatives[s[0]]['sequences']:
				if prev_seq not in key_sequences:
					key_sequences.append(prev_seq)
				sequences[prev_seq]['cluster_id'] = str(cluster)
				sequences[prev_seq]['identity_to_rep'] = s[1]	
			#updating sequences for previous representatives in sequences file
			sequences[s[0]]['cluster_id'] = str(cluster)
			sequences[s[0]]['identity_to_rep'] = s[1]
			for p in prev_key_sequences:
				sequences[p]['identity_to_rep'] = s[1]
				sequences[p]['cluster_id'] = str(cluster)
		cluster += 1  
		new_representatives[representative]['sequences'] = key_sequences

				
	#Getting representative sequences of previous step that are not included by blast in allvsall comparison
	pr = set(total_sequences)
	diff = [x for x in prev_representatives.keys() if x not in pr]
	print ( str(len(diff)) + '  clusters from previous step' )

	count = cluster
	for x in diff:
		new_representatives[x]= {}
		new_representatives[x]["cluster_id"] = str(count)
		#update cluster info in representative
		sequences[x]['cluster_id'] = str(count)
		sequences[x]['identity_to_rep'] = '100.00'
		new_representatives[x]['sequences'] = prev_representatives[x]['sequences']

		#update sequences in previous representatives
		for prev_seq in prev_representatives[x]['sequences']:
			sequences[prev_seq]['cluster_id'] = str(count)
		count=count+1

	print (str(count) + " clusters total identified for 70%")

	saveJson(new_representatives, 'dics/0.7repsgene.json')
	saveJson(sequences, 'dics/0.7seqsgene.json')

	
	#writing representatives file
	fasta_sequences = []
	print ("Writing representatives File")
	seqs_file = open(file_reps,'r')
	records = list(SeqIO.parse(seqs_file,'fasta'))
	for i in records:
		if (i.id in new_representatives.keys()): 
			fasta_sequences.append(i)
	output_handle = open('clustering/0.7-identity%s' %region_type, "w")
	SeqIO.write(fasta_sequences, output_handle, "fasta")		
	output_handle.close()

	buildMatrix(sequences, new_representatives, 0.7, samplesnumber, region_type)
	


if __name__ == "__main__":
	representatives_file = sys.argv[1]
	sequences_file = sys.argv[2]
	fasta_file = sys.argv[3]
	region_type = sys.argv[4]
	samples_file = sys.argv[5]
	allvsallclustering(representatives_file, sequences_file, fasta_file,region_type, samples_file)
	#allvsallclustering('../dics/0.95repsgene.json', '../dics/0.95seqsgene.json', '../clustering/0.95-identitygene','gene', '../dics/samples.json')