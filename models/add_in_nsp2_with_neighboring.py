import os
import collections
import pandas
import matplotlib, seaborn, numpy
from matplotlib import pyplot

import sklearn
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import SGDClassifier, LogisticRegression 
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score

import nsp2_interface
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


from sklearn.model_selection import train_test_split

import sklearn
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import SGDClassifier, LogisticRegression, LogisticRegressionCV
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score

import parse_bepipred

def count_aa(string):
	total_counts = []
	amino_acids = ['H', 'E', 'V', 'A', 'N', 'M', 'K', 'F', 'I', 'P', 'D', 'R', 'Y', 'T', 'S', 'W', 'G', 'C', 'L', 'Q']
	for aa in amino_acids:
		total_counts.append(string.lower().count(aa.lower()))
	return total_counts


nsp2 = nsp2_interface.NSP2()


def get_vectorized_sequences(list_sequences):
	vectorized_sequences = []
	for num in range(len(list_sequences)):
		sequence = list_sequences[num]
		num_aa = count_aa(sequence)
		#normalized_num_aa = [c/len(sequence) for c in num_aa]
		#num_hydrophobic = []
		#final_vector =  normalized_num_aa# + [count_attribute(sequence, "charged")] + [count_attribute(sequence, "polar")] + [count_attribute(sequence, "nonpolar")]
		#final_vector = [count_attribute(sequence, "charged")] + [count_attribute(sequence, "polar")] + [count_attribute(sequence, "nonpolar")]
		vectorized_sequences.append(num_aa)
	return vectorized_sequences


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
#for num in range(20000):
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
fraction_first_over_second_threshold = 0.75
uniprot_groups_to_use = []
for uniprot_group in list_of_uniprot_groups:
	total_hits = [seq.num_hits for seq in uniprot_group.list_of_sequences]
	position = [seq.start for seq in uniprot_group.list_of_sequences]
	#print(position)
	print(total_hits)
	if sorted(total_hits)[-1] > meaningful_threshold:
		top_uniprot_hits.append(sorted(total_hits)[-1])

		if len(total_hits) >= 3:# and sorted(total_hits)[-1] > meaningful_threshold:
			if total_hits.index(max(total_hits)) == 0:
				adjacent = total_hits[1]
			elif total_hits.index(max(total_hits)) == len(total_hits)-1:
				adjacent = total_hits[len(total_hits)-2]
			else:
				adjacent = sorted([total_hits[total_hits.index(max(total_hits))-1], total_hits[total_hits.index(max(total_hits))+1]])[1]
			fraction = sorted(total_hits)[-1] / ((adjacent + sorted(total_hits)[-1]))
			#print(fraction)
			#print(fraction <= fraction_first_over_second_threshold)
			fracion_first_over_second.append(fraction)
			if fraction <= fraction_first_over_second_threshold:
				#print("yayayayayay")
				uniprot_groups_to_use.append(uniprot_group)

list_of_binary_hits = []
just_amino_acids = []
just_nsp2_features = []
combined = []
print(len(uniprot_groups_to_use))
print("_____________________")
for uniprot_group in uniprot_groups_to_use:
	total_hits = [seq.num_hits for seq in uniprot_group.list_of_sequences]
	#print(total_hits)
	second_highest = sorted(total_hits)[-1]
	for sequence in uniprot_group.list_of_sequences:
		num_aa = count_aa(sequence.sequence)
		just_amino_acids.append(num_aa)
		nsp2_features = nsp2.get_summed_attributes_vector(sequence.uniprot_id, int(sequence.start), int(sequence.end))
		just_nsp2_features.append(nsp2_features)
		combined.append(num_aa + nsp2_features)
		if sequence.num_hits >= second_highest:
			list_of_binary_hits.append(True)
		else:
			list_of_binary_hits.append(False)




cross_val_score_args = {'cv': 5, 'n_jobs': -1, "scoring": "roc_auc", "groups": df.uniprot_accession.values}
scores = collections.OrderedDict()





score = cross_val_score(
    LogisticRegression(),
    just_amino_acids,
    df.binary_label, #y
    **cross_val_score_args)
scores["LogisticRegression kmers=1"] = (score)#, pipeline)
print("Just Amino Acids")
print(score)
print(numpy.mean(score))
print("")



score = cross_val_score(
    LogisticRegression(),
    just_nsp2_features,
    df.binary_label, #y
    **cross_val_score_args)
scores["LogisticRegression kmers=1"] = (score)#, pipeline)
print("Just NSP2")
print(score)
print(numpy.mean(score))
print("")

score = cross_val_score(
    LogisticRegression(),
    combined,
    df.binary_label, #y
    **cross_val_score_args)
scores["LogisticRegression kmers=1"] = (score)#, pipeline)
print("Combined")
print(score)
print(numpy.mean(score))



#sequences, classification = parse_bepipred.get_bepipred_data()


"""

for num in range(len(df["sequence_aa"])):
	print(df["sequence_aa"][num])
	print(df["binary_label"][num])
"""
#for seq in df["sequence_aa"]:
#	print(seq)