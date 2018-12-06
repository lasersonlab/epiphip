
import os
import collections
import pandas
import pandas as pd
import matplotlib, seaborn, numpy
from matplotlib import pyplot

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

class Sequence():
	def __init__(self, uniprot_id, start, end, sequence, num_hits):
		self.uniprot_id = uniprot_id
		self.start = start
		self.end = end
		self.sequence = sequence
		self.num_hits = num_hits

class UniprotGroup():
	def __init__(self, list_of_sequences):
		self.list_of_sequences = list_of_sequences
		sorted_seq = sorted(self.list_of_sequences, key=lambda sequence: sequence.start)
		#print(sorted_seq)
		i = 0
		aa_sequence = ""
		while i<len(list_of_sequences):
			aa_sequence = aa_sequence + list_of_sequences[i].sequence
			i = i+2
		if (len(list_of_sequences) % 2 == 0):
			index = int(list_of_sequences[-1].start) - int(list_of_sequences[-2].end)
			aa_sequence = aa_sequence + list_of_sequences[-1].sequence[index:]
		self.aa_sequence = aa_sequence

list_of_uniprot_ids = []

list_of_sequences = []
df = pandas.read_csv("{}/data/hits.binarized.fdr_0.15.w_metadata.csv".format(os.path.dirname(os.path.realpath(__file__))))
print(len(df.uniprot_accession))

for num in range(len(df.uniprot_accession)-4):
#for num in range(2000):
	print(num)
	#print(df.iloc[num])
	uniprot_id = df.uniprot_accession[num]
	sequence = df.sequence_aa[num]
	start = df.peptide_position[num].split("-")[0]
	end = df.peptide_position[num].split("-")[1]
	list_of_uniprot_ids.append(uniprot_id)
	num_hits = 0
	for number in df.iloc[num][4:]:
		num_hits = num_hits + int(number)
	list_of_sequences.append(Sequence(uniprot_id, start, end, sequence, num_hits))
list_of_uniprot_ids = list(set(list_of_uniprot_ids))


list_of_uniprot_groups = []
for uniprot_id in list_of_uniprot_ids:
	new_list_of_sequences = []
	for seq in list_of_sequences:
		if seq.uniprot_id == uniprot_id:
			new_list_of_sequences.append(seq)
	list_of_uniprot_groups.append(UniprotGroup(new_list_of_sequences))



summary_data = pd.DataFrame()
list_of_rows = []
for sequence in list_of_sequences:
	row = [sequence.uniprot_id, sequence.start, sequence.end, sequence.sequence, sequence.num_hits]
	list_of_rows.append(row)



df = pd.DataFrame(list_of_rows,columns=['uniprot_id','start', 'end', 'sequence', 'num_hits'])
df.to_csv("{}/data/summary_data.csv".format(os.path.dirname(os.path.realpath(__file__))))

meaningful_threshold = 10

non_zero_hits = [num for num in df.iloc[:,-1] if num > meaningful_threshold]
ax = sns.distplot(non_zero_hits,kde=False, rug=True)
fig = ax.get_figure()
fig.savefig("{}/data/non_zero_hits_per_sequence.png".format(os.path.dirname(os.path.realpath(__file__))))
plt.clf()

top_uniprot_hits = []
fracion_first_over_second = []
for uniprot_group in list_of_uniprot_groups:
	#print(uniprot_group.list_of_sequences[0].uniprot_id)
	#sorted_seq = sorted(uniprot_id.list_of_sequences, key=lambda sequence: sequence.num_hits)
	total_hits = [seq.num_hits for seq in uniprot_group.list_of_sequences]
	print(total_hits)
	print(sorted(total_hits)[-1])
	if sorted(total_hits)[-1] > meaningful_threshold:
		top_uniprot_hits.append(sorted(total_hits)[-1])
		if len(total_hits) > 1:# and sorted(total_hits)[-1] > meaningful_threshold:
			fracion_first_over_second.append(sorted(total_hits)[-1] / ((sorted(total_hits)[-2] + sorted(total_hits)[-1])))
		#print(sorted_seq)
		#print("{},{},{},{}".format(uniprot_id.list_of_sequences[0], uniprot_id.list_of_sequences[1]))#, uniprot_id.list_of_sequences[2], uniprot_id.list_of_sequences[3]))



ax = sns.distplot(top_uniprot_hits,kde=False, rug=True)
fig = ax.get_figure()
fig.savefig("{}/data/top_uniprot_hits.png".format(os.path.dirname(os.path.realpath(__file__))))
plt.clf()

ax = sns.distplot(fracion_first_over_second,kde=False, rug=True)
fig = ax.get_figure()
fig.savefig("{}/data/fracion_first_over_second.png".format(os.path.dirname(os.path.realpath(__file__))))
plt.clf()