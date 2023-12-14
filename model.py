from flaml import AutoML
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, recall_score, precision_score, confusion_matrix
import warnings
import preprocessing
import seaborn as sns
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


cols_to_encode = ['REF', 'ALT', 'GENE_SYMBOL', 'TYPE', 'VARIANT_TYPE', 
                      'MOST_SEVERE_CONSEQUENCE', 'IMPACT', 'EXON_INTRON_TYPE', 'CLINPRED_PRED']
cols_to_impute = ['CLINPRED_RANKSCORE', 'POLYPHEN2_HDIV_RANKSCORE',  'SIFT_CONVERTED_RANKSCORE',
                'SIFT4G_CONVERTED_RANKSCORE',  'MUTATIONASSESSOR_RANKSCORE', 'MUTATIONTASTER_CONVERTED_RANKSCORE']
cols_to_scale = ['STRAND', 'VARIANT_OCCURRENCES', 'EXON_INTRON_N', 'DOMAINS_COUNT', 'PUBMED_COUNT']


def load_data(path):

    df = pd.read_csv(path)
    df = df[df["LABEL"] != 2]
    y = df["LABEL"]
    x = df.drop(["LABEL"], axis=1)
    
    x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=0.2, stratify=y, random_state=47)
    return x_train, x_test, y_train, y_test



path = "data/dataset.csv"
x_train, x_test, y_train, y_test = load_data(path)

model = AutoML()
pipeline = preprocessing.Create_Pipeline(cols_to_encode, cols_to_scale, cols_to_impute=cols_to_impute)

x_train = pipeline.fit_transform(x_train)
model.fit(x_train, y_train, task="classification", time_budget=-1, verbose=0, metric="macro_f1", max_iter=100, eval_method="cv")

x_test = pipeline.transform(x_test)
pred = model.predict(x_test)

f_score = f1_score(y_test, pred, average="macro")
recall = recall_score(y_test, pred, average="macro")
precision = precision_score(y_test, pred, average="macro")

print(f"F1 score: {round(f_score, 4)} \nRecall score: {round(recall, 4)} \nPrecision: {round(precision, 4)} \nModel: {model.best_estimator}")
print()
confusion_matrix = confusion_matrix(y_test, pred) / len(y_test) * 100
sns.heatmap(confusion_matrix, annot=True)
plt.show()