# leak_demo.py -- 把"记忆"和"预测"分开：同一个数据集，同一个模型，只换划分方式。
# 模型是 1-NN（Tanimoto 最近邻），一个会背答案的模型。
# 关键不是它多准，而是：随机划分下的分数里，有多少是认脸认出来的。
import sys, json
sys.path.insert(0, ".")
import chemio, gates

rows = gates.load_csv("data/egfr_chembl203_ic50.csv")

def knn(tr, te):
    ts = [(r, r["bits"]) for r in tr]
    preds, ys, sims = [], [], []
    for r in te:
        best, bl = -1.0, None
        for rr, b in ts:
            t = chemio.tanimoto(r["bits"], b)
            if t > best: best, bl = t, rr
        preds.append(1.0 if bl["p"] >= 6.0 else 0.0)
        ys.append(1.0 if r["p"] >= 6.0 else 0.0)
        sims.append(best)
    n1 = sum(ys); n0 = len(ys) - n1
    pairs = sum(1 for i in range(len(ys)) for j in range(len(ys)) if ys[i]==1 and ys[j]==0 and preds[i]>preds[j])
    ties  = sum(1 for i in range(len(ys)) for j in range(len(ys)) if ys[i]==1 and ys[j]==0 and preds[i]==preds[j])
    auc = (pairs + 0.5*ties)/(n1*n0) if n1 and n0 else 0.5
    return auc, sum(sims)/len(sims)

out = {"model": "1-NN, Tanimoto on 2048-bit Morgan (a model that can memorise)"}
for name, (tr, te) in {
    "row_random":      gates.split_rows(rows),
    "unseen_molecule": gates.split_group(rows, lambda r: r["id"]),
    "unseen_document": gates.split_group(rows, lambda r: r["doc"]),
}.items():
    a, s = knn(tr, te)
    out[name] = {"auc": round(a,3), "mean_top1_tanimoto": round(s,3)}
out["gap_random_minus_unseen_molecule"] = round(out["row_random"]["auc"] - out["unseen_molecule"]["auc"], 3)
out["verdict"] = "REJECT" if out["gap_random_minus_unseen_molecule"] >= 0.10 else "PASS"
out["reason"] = ("随机划分下的分数来自认脸，不是来自化学"
                 if out["verdict"] == "REJECT" else "迁移：换分子/换论文后分数基本不掉")
print(json.dumps(out, ensure_ascii=False, indent=1))
