import os
import re
import sys
import subprocess
import re
import csv
import glob
import json
import itertools
import numpy as np
from pathlib import Path
from Bio import SeqIO

class Cluster:
	def __init__(self, cluster_id, representative, sequences):
		self.cluster_id = cluster_id
		self.representative = representative
		self.sequences = sequences

class Sequence:
	def __init__(self, sequence_id, cluster_id, sample_id, sequence_type):
		self.sequence_id = sequence_id
		self.cluster_id = cluster_id
		self.sample_id = sample_id
		self.distance_to_representative_cd_hit = 0.0
		self.distance_to_representative_mapping_blast = 0.0
		self.distance_to_representative_mcl_step = 0.0
		self.processed_cd_hit = False
		self.processed_mapping_blast = False
		self.processed_mcl = False
		self.sequence_type = sequence_type

def IterativeClusterSeqs(fasta_file, T, region_type, samples_file, results_name, global_file):

	#load samples File	
	samples_dic = readJson(samples_file)
	#len of samples
	samples = len(samples_dic.keys())
	
	sequences = {}
	clusters = {}
	representatives = {}
	percentages = [1, 0.99, 0.95, 0.9, 0.85, 0.8]
	prefix = '-identity%s' %region_type

	if region_type=='gene':
		n = 9
	else:
		n = 7

	iterative_step = False

	#Iterative Clustering
	for i in range(len(percentages)):

		identity = percentages[i] * 100
		mapped_blast = 0
		
		if i != 0:
			iterative_step = True
		
		#Executing CD-HIT program
		subprocess.call(['cd-hit-est', '-i', fasta_file, '-c', str(percentages[i]), '-o', 
			'clustering/%s' % str(percentages[i]) +prefix ,'-T', str(4), '-aL',str(percentages[i]), '-n', 
			str(n), '-d', '150','-r','1', '-M', '0'])
		
		#CD-HIT clustering Output File
		clustering_file = 'clustering/%s' % str(percentages[i]) +prefix+".clstr"
		generator_of_clustering_groups_step = []
			
		cluster_file = open(clustering_file)

		generator_of_clustering_groups_step = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))

		sequences = {}
		clusters = {}
		if iterative_step:
			prev_sequences = readJson('dics/%sseqs%s.json' % (str(percentages[i-1]),region_type))
			prev_clusters = readJson('dics/%sclstrs%s.json' % (str(percentages[i-1]),region_type))

		sequences_fasta = []  

		for cluster_result in generator_of_clustering_groups_step:

			#Cluster Name EX: >Cluster 1002
			cluster_name = cluster_result.__next__().strip()
			sequence_tuples = [(seq.split('>')[1].split('...')[0].split('|')[0],seq.split('>')[1].split('...')[0],seq.split('>')[1].split('...')[1]) for seq in generator_of_clustering_groups_step.__next__()]
			cluster_number = cluster_name.split(" ")[1]
			tuple_sequences = []
			representative = next(s[1] for s in sequence_tuples if '*' in s[2]) 

			representatives_number = 0
			sequences_number = 0
				
			for s in sequence_tuples:
				tuple_sequences.append(s[1])
				if s[1] != representative:
					sequences_fasta.append(s[1])
				if '*' in s[2]:
					x='100.00'
				if '*' not in s[2]:
					m = re.match('\s(at)\s(\S\S)(\d+)\.(\d+)',s[2])
					x=m.group(3)+'.'+m.group(4)
				sequence_object = Sequence(s[1], cluster_number, s[0], region_type)
				sequence_object.processed_cd_hit = True
				sequence_object.distance_to_representative_cd_hit = x
				sequences[s[1]] = sequence_object.__dict__
				if iterative_step:
					representative_cluster = next(item for item in prev_clusters.values() if item['representative'] == s[1])
					if representative_cluster:
						for seq in representative_cluster["sequences"]:
							#print ("Updating Sequences of cluster %s" %str(cluster_number) +" and prev rep %s" %s[1])
							sequence = prev_sequences[seq]
							sequence["cluster_id"] = cluster_number
							sequence["distance_to_representative_cd_hit"] = x
							sequences[sequence["sequence_id"]] = sequence
							if sequence["sequence_id"] not in tuple_sequences:
								tuple_sequences.append(sequence["sequence_id"])
							#If sequence is different of current representative and sequence is not in sequences fasta append it
							if sequence["sequence_id"] != representative and sequence["sequence_id"] not in sequences_fasta:
								sequences_fasta.append(sequence["sequence_id"])
					else:
						pass

					
			cluster_object = Cluster(cluster_number,representative, tuple_sequences)
			clusters[cluster_number] = cluster_object.__dict__

		if iterative_step:
			total_file = global_file
			seqs_file = open(total_file,'r')
			records = list(SeqIO.parse(seqs_file,'fasta'))

			representatives_step_file = 'clustering/reps%s.fasta' %str(percentages[i])
			representatives_ids = []
			output_handle1 = open(representatives_step_file, "w")
			print ("writing representatives file for Mapping with BLASTn")
			for c in clusters:
				rep_step = clusters[c]["representative"]
				seqs_step = clusters[c]["sequences"]
				if len(seqs_step) > 1:
					record_rep = next((item for item in records if item.id == rep_step), None)
					if record_rep:
						representatives_ids.append(record_rep)
					else:
						pass
			SeqIO.write(representatives_ids, output_handle1, "fasta")		
			output_handle1.close()

			sequences_step_file = 'clustering/seqs%s.fasta' %str(percentages[i])
			output_handle2 = open(sequences_step_file, "w")
			sequences_ids = []
			print ("writing sequences file for Mapping with BLASTn")
			for sequence_id in sequences_fasta:
				record_sequence = next((item for item in records if item.id == sequence_id),None)
				if record_sequence:
					sequences_ids.append(record_sequence)
				else:
					pass
			SeqIO.write(sequences_ids, output_handle2, "fasta")		
			output_handle2.close()

			representatives_number = len(representatives_ids)
			sequences_number = len(sequences_ids)

			output_step_file = 'clustering/blasted%s.fasta' %str(percentages[i])

			subprocess.call(['blastn', '-subject', representatives_step_file,'-query', sequences_step_file, '-out', output_step_file, '-num_alignments', '1', '-qcov_hsp_perc','100','-outfmt', "6 sseqid qseqid bitscore qlen length pident qcovs"])

			with open(output_step_file, 'r') as csvfile:
				comparisons = 0
				print ("reading blasted file")
				blastedreader = csv.reader(csvfile, delimiter='\t')
				for row in blastedreader:
					sequence_blast_step = sequences[row[1]]
					sequence_blast_id = sequence_blast_step['sequence_id']
					cluster = sequence_blast_step['cluster_id']
					representative_blast_step = clusters[cluster]['representative']
					comparisons += 1
					if row[1] == sequence_blast_id and row[0] == representative_blast_step:
						mapped_blast += 1
						sequences[sequence_blast_id]['distance_to_representative_mapping_blast'] = float(row[5])
						sequences[sequence_blast_id]['processed_mapping_blast'] = True
				print ("Comparisons")
				print (comparisons)
		#saving overview file
		total_sequences = len(sequences.items())
		total_clusters = len(clusters.items())
		unique_clusters = sum( 1 for item in clusters if len(clusters[item]['sequences']) == 1)

		with open(results_name, 'a') as f:
			f.write(str(identity) + "\t" + str(total_sequences) + "\t"+ str(total_clusters) + "\t" + region_type + "\t" + str(mapped_blast) + "\t\t" + str(representatives_number) + "\t" + str(sequences_number) + "\t" + str(unique_clusters) + "\n")

		buildMatrix(sequences, clusters, str(percentages[i]), samples, region_type)
		fasta_file = 'clustering/%s' % str(percentages[i]) + prefix
		saveJson(sequences, 'dics/%sseqs%s.json' % (str(percentages[i]),region_type))
		saveJson(clusters, 'dics/%sclstrs%s.json' % (str(percentages[i]), region_type))

def buildMatrix(sequences, clusters, perc_iter, samplesnumber, region_type):
	try:
		os.stat("matrix")
	except:
		os.mkdir("matrix")
	
	bin_name = 'matrix/%s_BIN_' % perc_iter + '%s.npy' %region_type
	freq_name = 'matrix/%s_FREQ_' % perc_iter + '%s.npy' %region_type
	perc_name = 'matrix/%s_PERC_' % perc_iter + '%s.npy' %region_type

	mbinary = np.zeros((samplesnumber,len(clusters)), int)
	mfreqs = np.zeros((samplesnumber,len(clusters)), int)
	mperc = np.zeros((samplesnumber,len(clusters)), float)


	cluster_table = {}
	for key in clusters:
		cluster_number =  key
		cluster_table[cluster_number] = {} 

		for s in clusters[key]['sequences']:
			sample_id = s.split('|')[0]
			if not sample_id in cluster_table[cluster_number].keys():
				cluster_table[cluster_number][sample_id] = {}
				cluster_table[cluster_number][sample_id]['sequences'] = [] 
				cluster_table[cluster_number][sample_id]['percentages'] = []
			cluster_table[cluster_number][sample_id]['sequences'].append(s)
			if float(sequences[s]['distance_to_representative_mcl_step']) > 0:
				cluster_table[cluster_number][sample_id]['percentages'].append(float(sequences[s]['distance_to_representative_mcl_step']))
			elif float(sequences[s]['distance_to_representative_mapping_blast']) > 0:
				cluster_table[cluster_number][sample_id]['percentages'].append(float(sequences[s]['distance_to_representative_mapping_blast']))
			else:
				cluster_table[cluster_number][sample_id]['percentages'].append(float(sequences[s]['distance_to_representative_cd_hit']))

	for key, value in cluster_table.items():
		for sp, dic in value.items():
			mbinary[int(sp.split('-')[1]),int(key)] = int(1)
			mfreqs[int(sp.split('-')[1]),int(key)] = len(dic['sequences'])
			mperc[int(sp.split('-')[1]),int(key)] = sum(dic['percentages']) / float(len(dic['percentages']))

	np.save(bin_name, mbinary)
	np.save(freq_name, mbinary)
	np.save(perc_name, mperc)
	
	saveJson(cluster_table, 'dics/%s%scluster_table.json' % (perc_iter,region_type))			


def saveJson(dic, filename):
	with open(filename, 'w') as f:
		json.dump(dic, f, indent=4)
		
def readJson(filename):
	with open(filename, 'r') as f:
		data = json.load(f)
	return data

if __name__ == "__main__":
	fasta_file = sys.argv[1]
	region_type = sys.argv[2]
	samples_file = sys.argv[3]
	results_name = '%sresults.txt' % region_type
	with open(results_name, 'w') as f:
		f.write("Iden %" + "\t" + "Seqs" + "\t"+ "Clstrs" + "\t" + "type" + "\t" + "mp_blast" + "\t" + "mp_rep" + "\t" "mp_seqs" + "\t" + "Unique" + "\n")
		f.write("-----" + "\t" + "-----" + "\t"+ "-----" + "\t" + "-----" + "\t" + "-----" + "\t\t" + "-----" + "\t" + "-----" + "\t" + "-----" + "\n")	 
	IterativeClusterSeqs(fasta_file, str(4), region_type, samples_file, results_name, fasta_file)
