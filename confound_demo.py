# confound_demo.py -- 真实数据上的泄露：把"批次"当成特征喂进去。
# assay 身份（同一实验条件下测出来的一批化合物）会带来批次效应：
# 模型学会了"这批化合物平均活性高"，随机划分下看着很准，
# 换成"没见过的 assay"划分就塌。四关在真数据上照样会响。
import sys, json
sys.path.insert(0, ".")
import gates

rows = gates.load_csv("data/egfr_chembl203_ic50.csv")

def with_assay(tr, te, K=64):
    """指纹 + 描述符 + assay one-hot（泄露通道）。"""
    assays = sorted(set(r["assay"] for r in tr))
    idx = {a: i for i, a in enumerate(assays)}
    Rtr, Rte = gates.featurize(tr, te, K)
    def add(R, rs):
        out = []
        for v, r in zip(R, rs):
            hot = [0.0] * len(assays)
            j = idx.get(r["assay"])
            if j is not None: hot[j] = 1.0
            out.append(v + hot)
        return out
    return add(Rtr, tr), add(Rte, te)

def run(tr, te):
    Rtr, Rte = with_assay(tr, te)
    th = gates.fit(Rtr, [1.0 if r["p"] >= 6.0 else 0.0 for r in tr])
    return gates.auc(th, Rte, [1.0 if r["p"] >= 6.0 else 0.0 for r in te])

tr, te = gates.split_rows(rows);            a = run(tr, te)
tr3, te3 = gates.split_group(rows, lambda r: r["assay"]); b = run(tr3, te3)
tr4, te4 = gates.split_group(rows, lambda r: r["doc"]);    c = run(tr4, te4)
print(json.dumps({
 "claim": "指纹 + 描述符 + assay 身份 -> 活性",
 "row_random":      {"auc": round(a, 3), "note": "测试行与训练集共享 assay 的比例很高"},
 "unseen_assay":    {"auc": round(b, 3)},
 "unseen_document": {"auc": round(c, 3)},
 "gap_random_minus_unseen_assay": round(a - b, 3),
 "verdict": "REJECT" if (a - b) >= 0.10 else "PASS",
 "attribution": "分数主要来自 assay 批次身份，不是化学结构" if (a - b) >= 0.10 else "化学结构本身可迁移",
}, ensure_ascii=False, indent=1))
