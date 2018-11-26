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




df = pandas.read_csv("{}/data/train.20171126.csv".format(os.path.dirname(os.path.realpath(__file__))))
df["binary_label"] = df.ratio > 0.5


cross_val_score_args = {'cv': 5, 'n_jobs': -1, "scoring": "roc_auc", "groups": df.uniprot_accession.values}
scores = collections.OrderedDict()

#print(len(df["sequence_aa"]))





pipeline = Pipeline([
    ('vect', CountVectorizer(analyzer="char", ngram_range=(1, 1))),
    ('clf', LogisticRegression()),
])
pipeline.fit(
    df.sequence_aa.values,
    df.binary_label)


print(pipeline.predict(df.sequence_aa.values))
roc_auc = sklearn.metrics.roc_auc_score(df.binary_label, pipeline.predict(df.sequence_aa.values))
print(roc_auc)


score = cross_val_score(
    pipeline,
    df.sequence_aa.values,
    df.binary_label,
    **cross_val_score_args)
print(score)
import ipdb ; ipdb.set_trace()


"""
pipeline = Pipeline([
    ('vect', CountVectorizer(analyzer="char", ngram_range=(1, 1))),
    ('clf', LogisticRegression()),
])
score = cross_val_score(
    pipeline,
    df.sequence_aa.values,
    df.binary_label,
    **cross_val_score_args)
scores["LogisticRegression kmers=1"] = (score, pipeline)
score

"""
"""

vectorized_sequences = []
for num in range(len(df["sequence_aa"])):
	sequence = df["sequence_aa"][num]
	num_aa = count_aa(sequence)
	#num_hydrophobic = []
	final_vector =  num_aa# + [count_attribute(sequence, "charged")] + [count_attribute(sequence, "polar")] + [count_attribute(sequence, "nonpolar")]
	#final_vector = [count_attribute(sequence, "charged")] + [count_attribute(sequence, "polar")] + [count_attribute(sequence, "nonpolar")]
	vectorized_sequences.append(final_vector)

print(numpy.shape(vectorized_sequences))
print(numpy.shape(df.binary_label))

logist = LogisticRegression(n_jobs= -1,solver='lbfgs').fit(vectorized_sequences, df.binary_label)
print(logist.score(vectorized_sequences, df.binary_label))

i = 0
for thing in logist.predict(vectorized_sequences):
	if thing==True:
		i=i+1
print(i)

blue=0
for thing in df.binary_label:
	if thing==True:
		blue=blue+1
print(blue)
print(logist.predict_proba(vectorized_sequences))
#roc_auc = sklearn.metrics.roc_auc_score(df.binary_label, logist.predict_proba(vectorized_sequences))
#print(roc_auc)

"""



#roc_auc = sklearn.metrics.roc_auc_score([True, True, False, False], [True, True, False, False])
#print(roc_auc)
"""
score = cross_val_score(
    LogisticRegression(),
    vectorized_sequences,
    df.binary_label, #y
    **cross_val_score_args)
scores["LogisticRegression kmers=1"] = (score)#, pipeline)

print(score)
print(numpy.mean(score))
"""

"""

for num in range(len(df["sequence_aa"])):
	print(df["sequence_aa"][num])
	print(df["binary_label"][num])
"""
#for seq in df["sequence_aa"]:
#	print(seq)