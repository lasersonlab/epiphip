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


nsp2 = nsp2_interface.NSP2()

df = pandas.read_csv("{}/data/train.20171126.csv".format(os.path.dirname(os.path.realpath(__file__))))
df["binary_label"] = df.ratio > 0.5


cross_val_score_args = {'cv': 5, 'n_jobs': -1, "scoring": "roc_auc", "groups": df.uniprot_accession.values}
scores = collections.OrderedDict()

print(len(df["sequence_aa"]))
print(len(df["uniprot_accession"]))
print(list(df))

sequences = df["sequence_aa"]
just_amino_acids = []
just_nsp2_features = []
combined = []
for num in range(len(sequences)):
#for num in range(10):
	uniprot_id = df['id'][num].split("|")[1]
	start = int(df['id'][num].split("|")[2].split("-")[0])
	end = int(df['id'][num].split("|")[2].split("-")[1])
	num_aa = count_aa(df["sequence_aa"][num])
	just_amino_acids.append(num_aa)
	#print("{}, {}".format(start, end))
	nsp2_features = nsp2.get_summed_attributes_vector(uniprot_id, start, end)
	just_nsp2_features.append(nsp2_features)
	combined_vector = num_aa + nsp2_features
	#print(final_vector)
	combined.append(combined_vector)

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