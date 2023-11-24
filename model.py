from flaml import AutoML
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, accuracy_score

def load_data(path):

    df = pd.read_csv(path)
    print(df.shape)
    df = df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"], keep="last")
    print(df.shape)
    y = df["RIS."]
    x = df.drop(["RIS."], axis=1)
   
    X_train, X_test, y_train, y_test = train_test_split(x,y,test_size=0.2, stratify=y)
    return X_train, X_test, y_train, y_test


def model_train(X_train, y_train):
    model = AutoML()
    model.fit(X_train, y_train, task="classification", time_budget=10, metric="macro_f1")
    pred = model.predict(X_test)
    print(f1_score(y_test, pred, average="macro"))

path = "data/dataset.csv"
X_train, X_test, y_train, y_test = load_data(path)
model_train(X_train, y_train)