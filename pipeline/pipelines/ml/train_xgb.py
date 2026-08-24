import pandas as pd, xgboost as xgb
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import roc_auc_score

df = pd.read_csv("features.csv")
df = df[df["label"].notna()]                     # 라벨 있는 행만
FEAT = ["model","ligand_iptm","iptm","ptm","plddt","gpde","has_clash",
        "pocket_size","pocket_n_models","p2rank_dist"]
X = pd.get_dummies(df[FEAT])                      # model → one-hot, 나머진 숫자 그대로
y = df["label"].astype(int)
g = df["target"]                                 # 타겟 단위로 train/test 분리(누수 방지)

tr, te = next(GroupShuffleSplit(test_size=0.2, n_splits=1, random_state=0).split(X, y, g))
clf = xgb.XGBClassifier(n_estimators=300, max_depth=5, learning_rate=0.05)
clf.fit(X.iloc[tr], y.iloc[tr])

print("AUC:", roc_auc_score(y.iloc[te], clf.predict_proba(X.iloc[te])[:,1]))
for f, imp in sorted(zip(X.columns, clf.feature_importances_), key=lambda x: -x[1])[:10]:
    print(f, round(imp, 3))                       # 어떤 피처가 중요했나

# python train_xgb.py