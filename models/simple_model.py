import os
import collections
import pandas
import matplotlib, seaborn, numpy
from matplotlib import pyplot

import sklearn
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import SGDClassifier, LogisticRegression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score

def count_aa(string):
	total_counts = []
	amino_acids = ['H', 'E', 'V', 'A', 'N', 'M', 'K', 'F', 'I', 'P', 'D', 'R', 'Y', 'T', 'S', 'W', 'G', 'C', 'L', 'Q']
	for aa in amino_acids:
		total_counts.append(string.lower().count(aa.lower()))
	return total_counts

def count_attribute(string, attribute):
	charged = ['R', 'K', 'D', 'E']
	polar = ['Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'W']
	nonpolar = ['A', 'I', 'L', 'M', 'F', 'V', 'P', 'G']
	aa_list = None

	if attribute == "charged":
		aa_list = charged
	elif attribute == "polar":
		aa_list = polar
	elif attribute == "nonpolar":
		aa_list = nonpolar
	else:
		raise(ValueError())

	aa_sum = 0
	for aa in aa_list:
		aa_sum = aa_sum + string.lower().count(aa.lower())
	return(aa_sum)


print(count_attribute("AAAAAAA", "polar"))



df = pandas.read_csv("{}/data/train.20171126.csv".format(os.path.dirname(os.path.realpath(__file__))))
df["binary_label"] = df.ratio > 0.5

"""
new_column = []
for num in range(len(df["sequence_aa"])):
	if df["binary_label"][num]:
		new_column.append(df["sequence_aa"][num] + "1")
	else:
		new_column.append(df["sequence_aa"][num] + "0")
df["cheat_sequence"] = new_column
"""

cross_val_score_args = {'cv': 5, 'n_jobs': -1, "scoring": "roc_auc", "groups": df.uniprot_accession.values}
scores = collections.OrderedDict()

print(len(df["sequence_aa"]))

cross_val_score_args = {'cv': 5, 'n_jobs': -1, "scoring": "roc_auc", "groups": df.uniprot_accession.values}
scores = collections.OrderedDict()


vectorized_sequences = []
for num in range(len(df["sequence_aa"])):
	sequence = df["sequence_aa"][num]
	num_aa = count_aa(sequence)
	#num_hydrophobic = []
	final_vector =  num_aa + [count_attribute(sequence, "charged")] + [count_attribute(sequence, "polar")] + [count_attribute(sequence, "nonpolar")]
	#final_vector = [count_attribute(sequence, "charged")] + [count_attribute(sequence, "polar")] + [count_attribute(sequence, "nonpolar")]
	vectorized_sequences.append(final_vector)

#vectorized_sequences = [count_aa(seq) for seq in df.sequence_aa.values]

"""
pipeline = Pipeline([
    ('vect', CountVectorizer(analyzer="char", ngram_range=(1, 1))),
    ('clf', LogisticRegression()),
])
"""

score = cross_val_score(
    LogisticRegression(),
    vectorized_sequences,
    df.binary_label, #y
    **cross_val_score_args)
scores["LogisticRegression kmers=1"] = (score)#, pipeline)

print(score)

"""

for num in range(len(df["sequence_aa"])):
	print(df["sequence_aa"][num])
	print(df["binary_label"][num])
"""
#for seq in df["sequence_aa"]:
#	print(seq)