import xgboost as xgb
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

def load_data(path):

    dataframe = pd.read_csv(path)
    x = dataframe.drop(dataframe.columns[-1], axis=1)
    y = dataframe.iloc[:,-1]
    X_train, X_test, y_train, y_test = train_test_split(x,y,test_size=0.2, random_state=10)
    return X_train, X_test, y_train, y_test


path = ""
X_train, X_test, y_train, y_test = load_data(path)


def model_train(X_train, X_test, y_train, y_test):
    model = xgb.train(X_train,y_train)

    accuracy_score(y_test,y_pred=model.predict(X_test))
    classification_report(y_test,y_pred=model.predict(X_test))
    confusion_matrix(y_test,y_pred=model.predict(X_test))