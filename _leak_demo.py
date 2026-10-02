
import sys, json, random
sys.path.insert(0, ".")
import gates
rows = gates.load_csv("data/egfr_chembl203_ic50.csv")
# 把每条记录复制一份到训练集：测试集里 100% 的分子训练时都见过
tr = rows[:]; te = rows[:]
Rtr, Rte = gates.featurize(tr, te)
th = gates.fit(Rtr, [1.0 if r["p"]>=6.0 else 0.0 for r in tr])
a = gates.auc(th, Rte, [1.0 if r["p"]>=6.0 else 0.0 for r in te])
print(json.dumps({"leak_100pct_auc": round(a,3),
                  "honest_unseen_molecule_auc": 0.893,
                  "delta": round(a-0.893,3)}))
