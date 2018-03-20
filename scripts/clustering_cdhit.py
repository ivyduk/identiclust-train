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


def IterativeClusterSeqs(fasta_file, T, region_type, samples_file):

	samples_dic = readJson(samples_file)
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

	iterative = False

	for i in range(len(percentages)):
		
		print ('reading ' + fasta_file)
		# check of exist a previous dictionary
		previous_file = Path('clustering/%s' % str(percentages[i]) +prefix)

		if previous_file.is_file():
			iterative = True
		
		subprocess.call(['cd-hit-est', '-i', fasta_file, '-c', str(percentages[i]), '-o', 
			'clustering/%s' % str(percentages[i]) +prefix ,'-T', str(4), '-aL',str(percentages[i]), '-n', 
			str(n), '-d', '150','-r','1', '-M', '0'])

		clusterFile = 'clustering/%s' % str(percentages[i]) +prefix+".clstr"
		RepFile = 'clustering/%s' % str(percentages[i]) +prefix
		cluster_groups_step = []
			
		cluster_file = open(clusterFile)

		cluster_groups_step = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))

		if not iterative:
			sequences = {}
			representatives = {}
		else:
			sequences = readJson('dics/%sseqs%s.json' % (str(percentages[i-1]),region_type))
			representatives = {}
			prev_representatives = readJson('dics/%sreps%s.json' % (str(percentages[i-1]),region_type)) 

		for cluster_result in cluster_groups_step:

			name = cluster_result.__next__().strip()
			seqs = [(seq.split('>')[1].split('...')[0].split('|')[0],seq.split('>')[1].split('...')[0],seq.split('>')[1].split('...')[1]) for seq in cluster_groups_step.__next__()]
			cluster = name.split(" ")[1]

			if not iterative:
				key_sequences = []
				representative = next(s[1] for s in seqs if '*' in s[2]) 
				representatives[representative] = {}
				representatives[representative]["cluster_id"] = cluster
				
				for s in seqs:
					key_sequences.append(s[1])
					sequences[s[1]] = {}
					sequences[s[1]]['sample_id'] = s[0]
					sequences[s[1]]['cluster_id'] = cluster
					if '*' in s[2]:
						x='100.00'
						sequences[s[1]]['identity_to_rep'] = x
					if '*' not in s[2]:
						m = re.match('\s(at)\s(\S\S)(\d+)\.(\d+)',s[2])
						x=m.group(3)+'.'+m.group(4)
						sequences[s[1]]['identity_to_rep'] = x
				representatives[representative]['sequences'] = key_sequences
			else:
				representative = next(s[1] for s in seqs if '*' in s[2])
				prev_key_sequences = prev_representatives[representative]['sequences'] 
				key_sequences = prev_key_sequences 
				representatives[representative] = {}
				representatives[representative]["cluster_id"] = cluster
				
				for s in seqs:
					if s[1] not in key_sequences:
						key_sequences.append(s[1])
					sequences[s[1]]['cluster_id'] = cluster
					if '*' not in s[2]:
						m = re.match('\s(at)\s(\S\S)(\d+)\.(\d+)',s[2])
						x=m.group(3)+'.'+m.group(4)
						sequences[s[1]]['identity_to_rep'] = x
						for p in prev_key_sequences:
							sequences[p]['identity_to_rep'] = x
							sequences[p]['cluster_id'] = cluster 
				representatives[representative]['sequences'] = key_sequences

		buildMatrix(sequences, representatives, str(percentages[i]), samples, region_type)
		fasta_file = 'clustering/%s' % str(percentages[i]) + prefix
		saveJson(sequences, 'dics/%sseqs%s.json' % (str(percentages[i]),region_type))
		saveJson(representatives, 'dics/%sreps%s.json' % (str(percentages[i]), region_type))

def buildMatrix(sequences, representatives, perc_iter, samplesnumber, region_type):
	try:
		os.stat("matrix")
	except:
		os.mkdir("matrix")
	
	bin_name = 'matrix/%s_BIN_' % perc_iter + '%s.npy' %region_type
	freq_name = 'matrix/%s_FREQ_' % perc_iter + '%s.npy' %region_type
	perc_name = 'matrix/%s_PERC_' % perc_iter + '%s.npy' %region_type

	mbinary = np.zeros((samplesnumber,len(representatives.keys())), int)
	mfreqs = np.zeros((samplesnumber,len(representatives.keys())), int)
	mperc = np.zeros((samplesnumber,len(representatives.keys())), float)


	cluster_table = {}
	for key, value in representatives.items():
		cluster_number =  value['cluster_id']
		cluster_table[cluster_number] = {} 

		for s in representatives[key]['sequences']:
			sample_id = s.split('|')[0]
			if not sample_id in cluster_table[cluster_number].keys():
				cluster_table[cluster_number][sample_id] = {}
				cluster_table[cluster_number][sample_id]['sequences'] = [] 
				cluster_table[cluster_number][sample_id]['percentages'] = []
			cluster_table[cluster_number][sample_id]['sequences'].append(s)
			cluster_table[cluster_number][sample_id]['percentages'].append(float(sequences[s]['identity_to_rep']))

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
		json.dump(dic, f)
		
def readJson(filename):
	with open(filename, 'r') as f:
		data = json.load(f)
	return data

if __name__ == "__main__":
	fasta_file = sys.argv[1]
	region_type = sys.argv[2]
	samples_file = sys.argv[3]
	IterativeClusterSeqs(fasta_file, str(4), region_type, samples_file)

