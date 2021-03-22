# -*- coding: utf-8 -*-
"""
Pycharm Editor
Create by zhwei 
Python：3.7.0
"""

import pandas as pd
from pathlib import Path
import joblib
from Features import all_feature
import numpy as np

th = 0.5 #阈值
clf_feature_order = ["AAC","AAI","CT5","CTD","DPC","GAAC","GDPC","NT5","PAAC"]


def get_base_proba(test_features,feature_index):
    test_features = test_features.values
    base_feature = []
    for idx,clf in zip(feature_index,clf_feature_order):
        features = test_features[:,idx]
        model = joblib.load(f'./Models/BaseModel/{clf}.m')
        base_proba = model.predict_proba(features)[:,-1]
        base_feature.append(base_proba)
    return np.array(base_feature).T


def Model_pred(fastafile,timestamp,Result_PATH):
    test_full_features, feature_index = all_feature(fastafile)
    base_feature = get_base_proba(test_full_features,feature_index)
    meta_clf = joblib.load('./Models/MetaModel/META.m')
    result = meta_clf.predict_proba(base_feature)[:,-1]


    df = pd.DataFrame(list(zip(list(test_full_features.index), list(result))))
    df.columns = ['Name','Probability']
    df["Class"] = df["Probability"].apply(lambda x: "APP" if x >=th else "non-APP") #根据阈值划定
    resultfile = Path(Result_PATH).joinpath(str(timestamp) + '.csv')
    df.to_csv(resultfile, index=False, header=True)
    return resultfile


if __name__ == '__main__':
    pass