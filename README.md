# evidence-gate

给生成层加一个裁判。它不生产分子、不排靶点、不打分——它回答一个问题：

**这条主张的依据，在扰动下还站得住吗。**

EDA 那套骨架（search / prove / perturb / cover）原样保留，只是把"门"从网表等价性
换成对预测型主张的审查。oracle 不是湿实验——湿实验一次不可逆、昂贵、无法当门用。
门用的是可复现的扰动：换划分、置换特征、置换标签。

## 四关

    G0 阴性对照  置换标签重训。分数不塌到 0.5，说明整条管线是坏的，
                 此时任何高分都不算数（先证伪自己，再谈别人的分数）
    G1 复现门    同一流程再跑一遍，数字逐位相同
    G2 分组门    分子级 / 论文级 / 时间级划分：换一个维度掏空测试集还站得住吗
    G3 泄露门    测试行里有多少条与训练集共享同一个分子 —— 这条指标比 AUC 本身更重要
    G4 归因门    在哪个维度上塌，塌多少，逐条列出；没测过的显式列为 uncovered

判决只与"分数是怎么来的"有关，与分数高低无关。

## 真数据实跑（2026-10-02，两个公共靶点，数字可复现）

ChEMBL 公共活性数据，同一套模型（10 个描述符 + 64 位 Morgan 指纹 -> 逻辑回归，
全部标准库实现，见 chemio.py）。一次 13 秒。同一份数据连跑三次，数字逐位一致。

靶点 A  EGFR (CHEMBL203)，IC50 — 2144 条 / 1451 个分子 / 118 篇文献 / 206 个 assay

    关卡                     划分              AUC
    G1  复现门               row-random        0.884   其中 35.3% 的测试行，训练时见过同一分子
    G2  换分子               unseen-molecule   0.890
    G3  换论文               unseen-document   0.843
    G4  换年代               post-2000 留出    0.785
    G5  换批次               unseen-assay      0.854
    G0  阴性对照             标签置换          0.501   （塌了，所以管线可信）
    -> PASS：分数在换分子/换论文/换年代/换批次后都站得住

靶点 B  hERG (CHEMBL240)，IC50 — 2275 条 / 1750 个分子 / 244 篇文献 / 273 个 assay

    关卡                     划分              AUC
    G1  复现门               row-random        0.749   27.2% 的测试行与训练集共享分子
    G2  换分子               unseen-molecule   0.709
    G3  换论文               unseen-document   0.645   <- 掉出 0.70
    G4  换年代               post-2000 留出    0.453   <- 掉到抛硬币
    G5  换批次               unseen-assay      0.702
    G0  阴性对照             标签置换          0.511
    -> REJECT：换论文 0.645 < 0.70；换年代 0.453 < 0.65

同一个模型、同一个流程、同一份代码：一个通过，一个驳回。

hERG 这条在文献默认的随机划分下报 0.749，看着是个能用的模型。四关走完才知道，
它的分数有相当一部分来自"记住了这 244 篇论文里的化合物"——换成没见过的年份
直接掉到 0.453，比抛硬币还差。**排的是可信度，不是分数。** 这正是生成层自己
不会说、也不会知道的那句话。

### 这一步里被门抓住的，是我自己的 bug

第一版 gates.py 里，时间划分用 `shuffle` 取一半，而年份只有 True/False 两个 key，
于是"哪一半当测试集"由随机数决定，同一份数据两次运行 AUC 不同（0.464 / 0.571）。
指纹里的 `hash()` 也是同一个病：Python 对 str 每进程加盐，Morgan 指纹每次都变。
改成 hold_out_true 与 FNV-1a 稳定哈希之后，三次连跑逐位一致。**G1 复现门
抓的第一条主张，是本仓库自己的。** 如果它连自己都不查，它凭什么查别人的。

## 跑

    python gates.py                                   # 默认 EGFR
    python gates.py data/herg_chembl240_ic50.csv      # hERG
    python gate.py       ; python leak_demo.py ; python confound_demo.py

    leak_demo.py       1-NN 最近邻（会背答案的模型）；随机 0.809 / 未见分子 0.774
    confound_demo.py   指纹 + assay 批次身份当特征；随机 0.901 / 未见 assay 0.851
    gate.py            合成数据，自带一条故意泄露的主张：0.703 被驳回，0.589 通过

依赖：无。标准库跑完全部。results.json 存着上面两次真跑的完整数字。

## 为什么不用 RDKit / numpy

这台机器上子进程里 numpy 的 BLAS 会分配内存失败
（`OpenBLAS error: Memory allocation still failed after 10 retries`），
`numpy.random` 直接 `MemoryError`。所以 SMILES 解析、Morgan 指纹、逻辑回归、AUC
全部用标准库重写了。这不是清高，是被环境逼的，也让整个仓库能拷到任何机器上直接跑。

## 文件

    chemio.py            SMILES -> 原子/邻接，10 个描述符，2048 位 Morgan 指纹，Tanimoto
    gates.py             四关主程序（真数据，任意 ChEMBL CSV）
    gate.py              合成对照（含故意泄露的主张，确认门会响）
    leak_demo.py         会背答案的模型 + 四关
    confound_demo.py    批次身份当特征的泄露
    fetch_chembl.py      ChEMBL 活性抓取（纯 urllib）
    results.json         两次真跑的完整数字
    data/                EGFR 与 hERG 活性 CSV（含 SMILES 与来源文献）

## 这层为什么该存在

生成层已经拥挤：分子生成、靶点排序、虚拟筛选，人人都在做。验证层是空的：
没人回答"这个分数可不可信、依据是什么、换一个测试集还成不成立"。
这正是 EDA 里已经做过的事——覆盖率数字好看但接线错了。

这里做的是**裁判**，不是**加速器**。它不告诉你哪个分子好，它告诉你哪个结论可以先信。

MIT License.

## Related work by the same author

The same claim -- *a number is meaningless until it is shown to survive its own
verification* -- is made and measured in other domains:

- [autoforge](https://github.com/ZhangYangyi03/autoforge) -- a tool's fitness, until an oracle outside the tool agrees
- [agentic-eda](https://github.com/ZhangYangyi03/agentic-eda) -- a circuit's area, until equivalence to the reference netlist is proven
- [debt-verify](https://github.com/ZhangYangyi03/debt-verify) -- a debt clause decision, until it survives the published revision record
- [tool-market](https://github.com/ZhangYangyi03/tool-market) -- a tool's liveness, until the hash chain says which revision is live
- [agent-safety-bench](https://github.com/ZhangYangyi03/agent-safety-bench) -- a model's safety compliance, measured rather than assumed
