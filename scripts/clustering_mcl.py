from clustering_cdhit import readJson, saveJson
import json
import os
import sys
import subprocess
import csv
import operator
from Bio import SeqIO
from clustering_cdhit import buildMatrix



def allvsallclustering(clusters, sequences, file_reps, region_type, samples_file, results_name, global_file): 
	
	samples_dic = readJson(samples_file)
	samplesnumber = len(samples_dic.keys())
	mapped_blast = 0
	comparisons_blast = []

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
	subprocess.call(['blastn', '-db',  '%s/%s%s' % (directoryallvsall,file_name,region_type), '-query', file_reps, '-perc_identity', '%s' % str(70),'-out', '%s/%sallblasted' % (directoryallvsall,region_type), '-qcov_hsp_perc','100', '-outfmt', "6 qseqid sseqid bitscore qlen length pident qcovs" ])
	#dictionary with length of each sequence used in BLAST coomparison
	tuples_len = {}
	#dictionary that stores tuples as keys of two sequences and the value is the identity calculated
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
			if row[0] != row[1]:
				comparisons_blast.append(row[0])
				comparisons_blast.append(row[1])
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
	new_clusters = {}
	prev_clusters = readJson(clusters)
	
	#Update Representatives file
	for rep, results in clusters_dic.items():
		if rep not in total_sequences:
			total_sequences.append(rep)
		seqs = results
		representative = rep
		prev_cluster_number = sequences[representative]['cluster_id']
		#sequences stored in previous cluster of representatvive
		prev_key_sequences = prev_clusters[prev_cluster_number]['sequences'] 
		key_sequences = prev_key_sequences 
		new_clusters[str(cluster)] = {}
		new_clusters[str(cluster)]["cluster_id"] = str(cluster)
		new_clusters[str(cluster)]["representative"] = representative

		#update cluster info in representative
		sequences[rep]['cluster_id'] = str(cluster)
		sequences[rep]['distance_to_representative_mcl_step'] = '100.00'
		sequences[rep]['processed_mcl'] = True

		if not results:
			for sequence in prev_key_sequences:
				sequences[sequence]['cluster_id'] = str(cluster)

		for s in seqs:	
			if s[0] not in key_sequences:
				key_sequences.append(s[0])
			#update sequences in previous representatives
			cluster_number_seq = sequences[s[0]]['cluster_id']
			for prev_seq in prev_clusters[cluster_number_seq]['sequences']:
				if prev_seq not in key_sequences:
					key_sequences.append(prev_seq)
				sequences[prev_seq]['cluster_id'] = str(cluster)
				sequences[prev_seq]['distance_to_representative_mcl_step'] = s[1]
				sequences[prev_seq]['processed_mcl'] = True	
			#updating sequences for previous representatives in sequences file
			sequences[s[0]]['cluster_id'] = str(cluster)
			sequences[s[0]]['distance_to_representative_mcl_step'] = s[1]
			sequences[s[0]]['processed_mcl'] = True	

			for p in prev_key_sequences:
				sequences[p]['distance_to_representative_mcl_step'] = s[1]
				sequences[p]['cluster_id'] = str(cluster)
				sequences[p]['processed_mcl'] = True
		new_clusters[str(cluster)]['sequences'] = key_sequences
		cluster += 1 
	
	pr = set(total_sequences) #new_representatives
	diff = []
	sequence_to_remapping = []
	for x in prev_clusters.keys():
		#Getting representative sequences of previous step that are not included by blast in allvsall comparison
		if prev_clusters[x]['representative'] not in pr and prev_clusters[x]['representative'] not in comparisons_blast:
			print (prev_clusters[x]['representative'])
			diff.append(prev_clusters[x]['representative'])
		#Getting representative sequences of previous step to do re mapping
		if prev_clusters[x]['representative'] not in pr:
			sequence_to_remapping.append(prev_clusters[x]['representative'])
	print ( str(len(diff)) + '  clusters from previous step excluded for comparison' )
	print ( str(len(sequence_to_remapping)) + '  clusters from previous step assigned to new cluster' )


	count = cluster
	for x in diff:
		total_sequences.append(x)
		new_clusters[count]= {}
		new_clusters[count]["cluster_id"] = str(count)
		new_clusters[count]["representative"] = x
		cluster_number_prev_rep = sequences[x]['cluster_id']
		#update cluster info in representative
		sequences[x]['cluster_id'] = str(count)
		sequences[x]['processed_mcl'] = True
		new_clusters[count]['sequences'] = prev_clusters[cluster_number_prev_rep]['sequences']

		#update sequences in previous clusters
		for prev_seq in prev_clusters[cluster_number_prev_rep]['sequences']:
			sequences[prev_seq]['cluster_id'] = str(count)
			sequences[prev_seq]['processed_mcl'] = True
		count=count+1

	print (str(count) + " clusters total identified for 70%")

	saveJson(new_clusters, 'dics/0.7clstrs%s.json' %region_type)
	saveJson(sequences, 'dics/0.7seqs%s.json' %region_type)

	
	#writing representatives file
	fasta_sequences = []
	print ("Writing representatives File")
	seqs_file = open(global_file,'r')
	records = list(SeqIO.parse(seqs_file,'fasta'))
	for i in records:
		if (i.id in total_sequences): 
			fasta_sequences.append(i)
	output_handle = open('clustering/0.7-identity%s' %region_type, "w")
	SeqIO.write(fasta_sequences, output_handle, "fasta")		
	output_handle.close()

	sequences_fasta = []
	representatives_step_file = 'clustering/reps%s.fasta' %str(70.0)
	representatives_ids = []
	output_handle1 = open(representatives_step_file, "w")
	print ("writing representatives file for Mapping with BLASTn")

	for c in new_clusters:
		rep_step = new_clusters[c]["representative"]
		seqs_step = new_clusters[c]["sequences"]
		if len(seqs_step) > 1:
			try:
				record_rep = next(item for item in records if item.id == rep_step)
			except StopIteration as e:
				print(e)
				break
			representatives_ids.append(record_rep)
			sequences_fasta += [sequence for sequence in seqs_step if (sequence != rep_step and sequence not in sequences_fasta)]


	SeqIO.write(representatives_ids, output_handle1, "fasta")		
	output_handle1.close()

	sequences_step_file = 'clustering/seqs%s.fasta' %str(70.0)
	output_handle2 = open(sequences_step_file, "w")
	sequences_ids = []
	print ("writing sequences file for Mapping with BLASTn")
	for sequence_id in sequences_fasta:
		try:
			record_sequence = next(item for item in records if item.id == sequence_id)
			sequences_ids.append(record_sequence)
		except StopIteration as e:
			print(e)
			break
	SeqIO.write(sequences_ids, output_handle2, "fasta")		
	output_handle2.close()

	representatives_number = len(representatives_ids)
	sequences_number = len(sequences_ids)

	output_step_file = 'clustering/blasted%s.fasta' %str(70.0)

	subprocess.call(['blastn', '-subject', representatives_step_file,'-query', sequences_step_file, '-out', output_step_file, '-num_alignments', '1', '-qcov_hsp_perc','100','-outfmt', "6 sseqid qseqid bitscore qlen length pident qcovs"])

	with open(output_step_file, 'r') as csvfile:
		comparisons = 0
		print ("reading blasted file")
		blastedreader = csv.reader(csvfile, delimiter='\t')
		for row in blastedreader:
			sequence_blast_step = sequences[row[1]]
			sequence_blast_id = sequence_blast_step['sequence_id']
			cluster = sequence_blast_step['cluster_id']
			representative_blast_step = new_clusters[cluster]['representative']
			comparisons += 1
			if row[1] == sequence_blast_id and row[0] == representative_blast_step:
				mapped_blast += 1
				sequences[sequence_blast_id]['distance_to_representative_mapping_blast'] = float(row[5])
				sequences[sequence_blast_id]['processed_mapping_blast'] = True
		print ("Comparisons")
		print (comparisons)	

	len_total_sequences = len(sequences.items())
	total_clusters = len(total_sequences)
	unique_clusters = sum( 1 for item in new_clusters if len(new_clusters[item]['sequences']) == 1)



	with open(results_name, 'a') as f:
			f.write(str(70.0) + "\t" + str(len_total_sequences) + "\t"+ str(total_clusters) + "\t" + region_type + "\t" + str(mapped_blast) + "\t\t" + str(representatives_number) + "\t" + str(sequences_number) + "\t" + str(unique_clusters) + "\n")

	buildMatrix(sequences, new_clusters, 0.7, samplesnumber, region_type)
	


if __name__ == "__main__":
	clusters_file = sys.argv[1]
	sequences_file = sys.argv[2]
	fasta_file = sys.argv[3]
	region_type = sys.argv[4]
	samples_file = sys.argv[5]
	global_file = sys.argv[6]
	results_name = '%sresults.txt' % region_type
	allvsallclustering(clusters_file, sequences_file, fasta_file,region_type, samples_file, results_name, global_file)