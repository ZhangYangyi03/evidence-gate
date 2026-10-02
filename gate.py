# evidence-gate: 验证层最小闭环（eda-spine 四段骨架 -> "预测型主张"的裁判）。
# 批量模式：一次审 N 条主张，每条过四关，输出存活/驳回排名。
import random, json, math, time, sys

t0 = time.time()
G, PER, D = 24, 16, 10
random.seed(7)
g = [i for i in range(G) for _ in range(PER)]
N = G * PER
bias = [random.gauss(0, 2.2) for _ in range(G)]          # 骨架带来的强批次效应 = 泄露源
X = [[random.gauss(0, 1) for _ in range(D)] for _ in range(N)]
W = [0.6, 0.4, 0.3] + [0.0] * (D - 3)
def lab(x, b): return 1.0 if random.random() < 1/(1+math.exp(-(x+b))) else 0.0
y = [lab(sum(X[i][j]*W[j] for j in range(D)), bias[g[i]]) for i in range(N)]

def build(idxs, leaky, perm=None):
    R = []
    for i in idxs:
        v = list(X[i])
        if leaky:
            src = perm[i] if perm is not None else i
            v += [1.0 if g[src] == k else 0.0 for k in range(G)]
        R.append(v)
    return R

def fit(R, Y, iters=90, lr=0.7, l2=1e-3):
    m = len(R[0]) + 1; th = [0.0]*m; n = len(Y)
    for _ in range(iters):
        gr = [0.0]*m
        for r, yy in zip(R, Y):
            z = th[-1]
            for j in range(m-1): z += r[j]*th[j]
            d = 1/(1+math.exp(-max(-30.0, min(30.0, z)))) - yy
            for j in range(m-1): gr[j] += d*r[j]
            gr[-1] += d
        for j in range(m): th[j] -= lr*(gr[j]/n + l2*th[j])
    return th

def auc(th, R, Y):
    s = [th[-1] + sum(r[j]*th[j] for j in range(len(r))) for r in R]
    o = sorted(range(len(s)), key=lambda i: s[i]); rk = [0]*len(s)
    for p, i in enumerate(o): rk[i] = p+1
    n1 = sum(Y); n0 = len(Y)-n1
    if n1 == 0 or n0 == 0: return 0.5
    return (sum(rk[i] for i in range(len(Y)) if Y[i] == 1) - n1*(n1+1)/2)/(n1*n0)

def split(kind, seed=0, frac=0.5):
    r = random.Random(seed)
    if kind == "random":
        idx = list(range(N)); r.shuffle(idx); c = int(N*frac); return idx[:c], idx[c:]
    gs = list(range(G)); r.shuffle(gs); te = set(gs[:int(G*frac)])
    return [i for i in range(N) if g[i] not in te], [i for i in range(N) if g[i] in te]

def run(kind, leaky, perm=None, seed=0):
    tr, te = split(kind, seed)
    th = fit(build(tr, leaky, perm), [y[i] for i in tr])
    return auc(th, build(te, leaky, perm), [y[i] for i in te])

def audit(claim_name, leaky):
    rep = {"claim": claim_name, "gates": {}}
    a1, a2 = run("random", leaky), run("random", leaky)
    rep["gates"]["G1_reproduce"] = {"auc": round(a1, 3), "deterministic": abs(a1-a2) < 1e-12}
    base = a1
    drops = {}
    for j in range(D):
        perm = list(range(N)); random.Random(100+j).shuffle(perm)
        col = [X[i][j] for i in range(N)]
        Xs = [list(X[i]) for i in range(N)]
        for i in range(N): Xs[i][j] = col[perm[i]]
        Xsave = X[:]
        X[:] = Xs
        drops["desc%d" % j] = round(base - run("random", leaky), 3)
        X[:] = Xsave
    if leaky:
        perm = list(range(N)); random.Random(555).shuffle(perm)
        drops["scaffold_block"] = round(base - run("random", True, perm), 3)
    perm = list(range(N)); random.Random(777).shuffle(perm)
    ymix = [y[perm[i]] for i in range(N)]
    tr, te = split("random", 0)
    th = fit(build(tr, leaky), [ymix[i] for i in tr])
    negdrop = round(base - auc(th, build(te, leaky), [ymix[i] for i in te]), 3)
    top = sorted(drops.items(), key=lambda kv: -kv[1])[:2]
    rep["gates"]["G2_perturb"] = {"auc": round(base, 3), "top_drops": dict(top),
                                  "label_permute_drop": negdrop,
                                  "sanity_ok": negdrop >= 0.1}
    rnd, scf = run("random", leaky), run("scaffold", leaky)
    rep["gates"]["G3_split"] = {"random": round(rnd, 3), "scaffold": round(scf, 3),
                                "gap": round(rnd-scf, 3), "robust": (rnd-scf) < 0.10}
    space = ["split:random", "split:scaffold", "label_permute"] + ["permute:"+k for k in drops]
    ran = space[:] if leaky else space[:3] + ["permute:"+k for k in drops]
    rep["gates"]["G4_cover"] = {"space": len(space), "ran": len(ran),
                                "uncovered": [s for s in space if s not in ran]}
    g1 = rep["gates"]["G1_reproduce"]["deterministic"]
    g2 = rep["gates"]["G2_perturb"]["sanity_ok"]
    g3 = rep["gates"]["G3_split"]["robust"]
    rep["verdict"] = "PASS" if (g1 and g2 and g3) else "REJECT"
    rep["reason"] = "; ".join(x for x, ok in [("G1 不可复现", g1), ("G2 置换标签不塌(管线不可信)", g2),
                                              ("G3 换划分即塌(claim 不迁移)", g3)] if not ok) or "四关全过"
    return rep

claims = [("A: 描述符+骨架身份 -> 活性", True), ("B: 仅描述符 -> 活性", False)]
res = [audit(n, l) for n, l in claims]
print(json.dumps({"batch": res, "seconds": round(time.time()-t0, 1)}, ensure_ascii=False, indent=1))
