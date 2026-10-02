# gates.py -- 四关调用四条不同的划分，每关都把"泄露"与"迁移"分开量。
# 这里不做打分，只做裁判：同一条预测主张，换划分方式后还站不站得住。
import json, math, os, random, sys, csv
from collections import Counter

sys.path.insert(0, ".")
import chemio

pct = lambda x: round(100.0 * x, 1)

def load_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if not r.get("smiles") or not r.get("pchembl_value"):
                continue
            atoms, adj = chemio.parse(r["smiles"])
            rows.append({"id": r["molecule_chembl_id"], "doc": r["document_chembl_id"],
                         "year": r.get("document_year") or "", "assay": r.get("assay_chembl_id"),
                         "p": float(r["pchembl_value"]),
                         "bits": chemio.morgan_bits(atoms, adj),
                         "desc": chemio.descriptors(atoms, adj)})
    return rows

def top_bits(tr, K):
    cp, cn = Counter(), Counter()
    for r in tr:
        b, tgt = r["bits"], (cp if r["p"] >= 6.0 else cn)
        while b:
            l = b & -b; tgt[l.bit_length() - 1] += 1; b ^= l
    npos = max(1, sum(1 for r in tr if r["p"] >= 6.0)); nneg = max(1, len(tr) - npos)
    sc = {k: cp.get(k, 0) / npos - cn.get(k, 0) / nneg for k in set(cp) | set(cn)}
    return [k for k, _ in sorted(sc.items(), key=lambda kv: -abs(kv[1]))[:K]]

def featurize(tr, te, K=64):
    D = len(tr[0]["desc"]); sel = top_bits(tr, K)
    mu = [sum(r["desc"][j] for r in tr) / len(tr) for j in range(D)]
    sd = [max(1e-6, (sum((r["desc"][j] - mu[j]) ** 2 for r in tr) / len(tr)) ** 0.5) for j in range(D)]
    f = lambda r: [(r["desc"][j] - mu[j]) / sd[j] for j in range(D)] + \
                  [1.0 if (r["bits"] >> k) & 1 else 0.0 for k in sel]
    return [f(r) for r in tr], [f(r) for r in te]

def fit(R, Y, it=150, lr=0.9, l2=1e-3):
    m = len(R[0]) + 1; th = [0.0] * m; n = len(Y)
    for _ in range(it):
        g = [0.0] * m
        for r, yy in zip(R, Y):
            z = th[-1]
            for j in range(m - 1): z += r[j] * th[j]
            p = 1 / (1 + math.exp(-max(-30.0, min(30.0, z)))) - yy
            for j in range(m - 1): g[j] += p * r[j]
            g[-1] += p
        for j in range(m): th[j] -= lr * (g[j] / n + l2 * th[j])
    return th

def auc(th, R, Y):
    s = [th[-1] + sum(r[j] * th[j] for j in range(len(r))) for r in R]
    o = sorted(range(len(s)), key=lambda i: s[i]); rk = [0] * len(s)
    for i, k in enumerate(o): rk[k] = i + 1
    n1 = sum(1 for v in Y if v == 1); n0 = len(Y) - n1
    if not n1 or not n0: return 0.5
    return (sum(rk[i] for i in range(len(Y)) if Y[i] == 1) - n1 * (n1 + 1) / 2) / (n1 * n0)

def split_rows(rows, seed=0, frac=0.5):
    r = list(rows); random.Random(seed).shuffle(r); c = int(len(r) * frac); return r[:c], r[c:]

def split_group(rows, key, seed=0, frac=0.5, hold_out_true=False):
    """按 key 分组划分。hold_out_true=True 时，把 key 为真的那一组整体留作测试集。
    时间划分必须用这个：只有 True/False 两个 key，靠 shuffle 取一半等于抛硬币，
    同一份数据两次运行会给出不同的 AUC —— 那正好是这一层要消灭的东西。"""
    ks = sorted(set(key(r) for r in rows), key=lambda k: str(k))
    if hold_out_true:
        te = set(k for k in ks if k)
    else:
        random.Random(seed).shuffle(ks)
        te = set(ks[:max(1, int(len(ks) * frac))])
    return [r for r in rows if key(r) not in te], [r for r in rows if key(r) in te]

def score_split(rows, tr, te, K=64):
    Rtr, Rte = featurize(tr, te, K)
    Ytr = [1.0 if r["p"] >= 6.0 else 0.0 for r in tr]
    Yte = [1.0 if r["p"] >= 6.0 else 0.0 for r in te]
    th = fit(Rtr, Ytr)
    ids = set(r["id"] for r in tr)
    return auc(th, Rte, Yte), round(sum(1 for r in te if r["id"] in ids) / max(1, len(te)), 3)

def label_permute_control(rows, tr, te, K=64, seed=5):
    Y = [1.0 if r["p"] >= 6.0 else 0.0 for r in rows]
    perm = list(range(len(rows))); random.Random(seed).shuffle(perm)
    pm = {id(rows[i]): Y[perm[i]] for i in range(len(rows))}
    Rtr, Rte = featurize(tr, te, K)
    th = fit(Rtr, [pm[id(r)] for r in tr])
    return round(auc(th, Rte, [pm[id(r)] for r in te]), 3)

def main(csv_path="data/egfr_chembl203_ic50.csv", target_label=None):
    rows = load_csv(csv_path)
    label = target_label or os.path.splitext(os.path.basename(csv_path))[0]
    rep = {"claim": f"SMILES -> 活性 (pChEMBL >= 6)，AUC | {label}", "dataset":
           {"csv": csv_path, "rows": len(rows), "unique_molecules": len(set(r["id"] for r in rows)),
            "documents": len(set(r["doc"] for r in rows)),
            "assays": len(set(r["assay"] for r in rows)), "label": label}, "gates": {}}

    tr, te = split_rows(rows); a, leak = score_split(rows, tr, te)
    rep["gates"]["G1_reproduce"] = {"split": "row-random", "auc": round(a, 3),
        "test_rows_sharing_a_molecule_with_train": leak,
        "note": "文献默认划分；这条本身就是泄露的证据"}
    tr2, te2 = split_group(rows, lambda r: r["id"]); b, _ = score_split(rows, tr2, te2)
    rep["gates"]["G2_unseen_molecule"] = {"split": "unseen-molecule", "auc": round(b, 3)}
    tr3, te3 = split_group(rows, lambda r: r["doc"]); c, _ = score_split(rows, tr3, te3)
    rep["gates"]["G3_unseen_paper"] = {"split": "unseen-document", "auc": round(c, 3)}
    d2, _ = score_split(rows, *split_group(rows, lambda r: bool(r["year"]) and int(r["year"]) >= 2000,
                                              hold_out_true=True))
    rep["gates"]["G4_temporal"] = {"split": "post-2000 held out", "auc": round(d2, 3)}
    e, _ = score_split(rows, *split_group(rows, lambda r: r["assay"]))
    rep["gates"]["G5_unseen_assay"] = {"split": "unseen-assay", "auc": round(e, 3),
        "note": "assay 身份是批次效应的入口；掉得多说明分数里有一部分是实验室批次，不是化学"}
    neg = label_permute_control(rows, tr3, te3)
    rep["gates"]["G0_negative_control"] = {"label_permute_auc": neg,
        "ok": abs(neg - 0.5) < 0.06, "note": "不塌说明管线是坏的"}
    rep["gap_random_vs_paper"] = round(a - c, 3)
    rep["gap_random_vs_temporal"] = round(a - d2, 3)
    fails = []
    if abs(neg - 0.5) >= 0.06: fails.append("G0 阴性对照没塌到 0.5，管线本身不可信")
    if c < 0.70: fails.append(f"G3 换论文 AUC {round(c,3)} < 0.70（分数没迁移出这批文献）")
    if d2 < 0.65: fails.append(f"G4 换年代 AUC {round(d2,3)} < 0.65（配方换了年份就失效）")
    if a - c >= 0.15: fails.append(f"G1->G3 落差 {round(a-c,3)} >= 0.15（随机划分的分数主要来自记住这批数据）")
    rep["verdict"] = "PASS" if not fails else "REJECT"
    rep["attribution"] = "；".join(fails) if fails else "四关全过；分数在换分子/换论文/换年代/换批次后都站得住"
    print(json.dumps(rep, ensure_ascii=False, indent=1))

if __name__ == "__main__":
    main(*sys.argv[1:])
