import os
import collections
import pandas
import matplotlib, seaborn, numpy
from matplotlib import pyplot

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



df = pandas.read_csv("{}/data/train.20171126.csv".format(os.path.dirname(os.path.realpath(__file__))))
df["binary_label"] = df.ratio > 0.5



vectorized_sequences = get_vectorized_sequences(df["sequence_aa"])

model = LogisticRegression()
model.fit(vectorized_sequences, df.binary_label.values)

roc_auc = sklearn.metrics.roc_auc_score(df.binary_label.values, model.predict_proba(vectorized_sequences)[:,1])
print("Roc score on training data: {}".format(roc_auc))


sequences, classification = parse_bepipred.get_bepipred_data()
print(sequences[0])
print(classification[0])
bepi_vectorized = get_vectorized_sequences(sequences)
bepi_roc_auc = sklearn.metrics.roc_auc_score(classification, model.predict_proba(bepi_vectorized)[:,1])
print("Roc score on bepipred data: {}".format(bepi_roc_auc))


from sklearn.model_selection import train_test_split


combined_data = bepi_vectorized
combined_classification = classification

X_train, X_test, y_train, y_test = train_test_split(combined_data, combined_classification, test_size=0.9, random_state=0)

scores = cross_val_score(LogisticRegression(solver='lbfgs'), X_train, y_train, cv=5, scoring='roc_auc')
print("cross_val_score for bepipred: {}".format(scores))
model = LogisticRegression(solver='lbfgs')
model.fit(X_train, y_train)

roc_auc = sklearn.metrics.roc_auc_score(y_test, model.predict_proba(X_test)[:,1])
print("Roc score on test bepipred data alone: {}".format(roc_auc))


combined_data = vectorized_sequences
combined_classification = [a for a in df.binary_label.values]
X_train, X_test, y_train, y_test = train_test_split(combined_data, combined_classification, test_size=0.9, random_state=0)

scores = cross_val_score(LogisticRegression(solver='lbfgs'), X_train, y_train, cv=5, scoring='roc_auc')
print("cross_val_score for phip-seq: {}".format(scores))
model = LogisticRegression(solver='lbfgs')
model.fit(X_train, y_train)

roc_auc = sklearn.metrics.roc_auc_score(y_test, model.predict_proba(X_test)[:,1])
print("Roc score on test phip-seq data alone: {}".format(roc_auc))




combined_data = vectorized_sequences + bepi_vectorized
combined_classification = [a for a in df.binary_label.values] + classification

X_train, X_test, y_train, y_test = train_test_split(combined_data, combined_classification, test_size=0.9, random_state=0)

scores = cross_val_score(LogisticRegression(solver='lbfgs'), X_train, y_train, cv=5, scoring='roc_auc')
print("cross_val_score for combined: {}".format(scores))
model = LogisticRegression(solver='lbfgs')
model.fit(X_train, y_train)

roc_auc = sklearn.metrics.roc_auc_score(y_test, model.predict_proba(X_test)[:,1])
print("Roc score on test data for combined: {}".format(roc_auc))

