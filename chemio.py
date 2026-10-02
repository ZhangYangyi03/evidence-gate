# chemio.py -- 纯标准库 SMILES 解析 + 描述符 + Morgan 指纹 + Tanimoto。
# 故意不依赖 RDKit / numpy：本机 BLAS 在子进程里分配内存失败（OpenBLAS 10 retries），
# 所以任何依赖 numpy 的路线在这台机器上都跑不起来。
MASS = {"C":12.011,"N":14.007,"O":15.999,"S":32.06,"F":18.998,"Cl":35.45,"Br":79.904,
        "I":126.90,"P":30.974,"B":10.81,"Si":28.086,"Se":78.97}
AROM = set("bcnops")

def parse(smiles):
    """-> (atoms, adj)  atoms: [element]; adj: {i: {j: bond_order}}  抛异常表示解析不了。"""
    atoms=[]; adj={}; cur=None; stack=[]; pend=None; rings={}
    i=0; n=len(smiles)
    while i<n:
        ch=smiles[i]
        if ch=="[":
            j=smiles.index("]", i)
            inner=smiles[i+1:j]
            m="".join(c for c in inner if c.isalpha()) or "*"
            el=(m[0].upper()+m[1:].lower()) if m[0].isupper() else m[0].upper()
            arom = m[0].islower()
            atoms.append((el,arom)); k=len(atoms)-1
            if pend is None: b=1.5 if (arom and cur is not None and atoms[cur][1]) else 1
            else: b={"-":1,"=":2,"#":3,":":1.5,"/":1,"\\":1}.get(pend,1); pend=None
            adj.setdefault(k,{})
            if cur is not None: adj[cur][k]=b; adj[k][cur]=b
            cur=k; i=j+1; continue
        if ch.isalpha():
            if i+1<n and smiles[i+1].islower() and ch.isupper() and (ch+smiles[i+1]).lower() in ("cl","br","si","se"):
                el=ch+smiles[i+1]; arom=False; i+=2
            else:
                el=ch.upper(); arom=ch.islower(); i+=1
            atoms.append((el,arom)); k=len(atoms)-1
            if pend is None: b=1.5 if (arom and cur is not None and atoms[cur][1]) else 1
            else: b={"-":1,"=":2,"#":3,":":1.5,"/":1,"\\":1}.get(pend,1); pend=None
            adj.setdefault(k,{})
            if cur is not None: adj[cur][k]=b; adj[k][cur]=b
            cur=k; continue
        if ch in "-=#:/\\": pend=ch; i+=1; continue
        if ch=="(": stack.append(cur); i+=1; continue
        if ch==")": cur=stack.pop(); i+=1; continue
        if ch==".": cur=None; i+=1; continue
        if ch.isdigit() or ch=="%":
            if ch=="%": num=int(smiles[i+1:i+3]); i+=3
            else: num=int(ch); i+=1
            b = 1.5 if (cur is not None and atoms[cur][1]) else 1
            if pend: b={"-":1,"=":2,"#":3,":":1.5}.get(pend,1); pend=None
            if num in rings:
                o=rings.pop(num); adj.setdefault(cur,{})[o]=b; adj.setdefault(o,{})[cur]=b
            else: rings[num]=cur
            continue
        if ch in "@+": i+=1; continue
        i+=1
    if not atoms: raise ValueError("empty")
    return atoms, adj

def inchi_key_like(atoms, adj):
    """没有 RDKit，用一个稳定的规范序近似：度/元素/环基分层 BFS 后的哈希。"""
    h=[]
    for k in range(len(atoms)):
        nb=sorted(adj.get(k,{}).items())
        h.append((len(nb), atoms[k][0], atoms[k][1], tuple(b for _,b in nb)))
    return "@".join(map(str,h))

def descriptors(atoms, adj):
    """10 个可解释的描述符，全部整数/浮点，无依赖。"""
    nH=len(atoms); deg={k:len(adj.get(k,{})) for k in range(nH)}
    rings = len(adj and [1]) and (sum(deg.values())//2 - nH + 1)   # 键-原子+1（单片段）
    arom=sum(1 for e,a in atoms if a)
    het=sum(1 for e,_ in atoms if e in ("N","O","S","P","F","Cl","Br","I"))
    hal=sum(1 for e,_ in atoms if e in ("F","Cl","Br","I"))
    mw=sum(MASS.get(e,12.0) for e,_ in atoms)
    rot=0
    for k in range(nH):
        for j,b in adj.get(k,{}).items():
            if j>k and b==1 and deg[k]>1 and deg[j]>1 and rings<=0: rot+=1
    tpsa=sum({"N":12.0,"O":17.0,"S":25.0,"P":23.0}.get(e,0.0) for e,_ in atoms)
    hbd=sum(1 for e,_ in atoms if e in ("N","O"))
    sp3=sum(1 for k in range(nH) if deg[k]<=3 and not atoms[k][1])/max(1,nH)
    return [nH, rings, arom, het, hal, round(mw,1), rot, round(tpsa,1), hbd, round(sp3,3)]

def morgan_bits(atoms, adj, radius=2, nbits=2048):
    """ECFP 风格：原子不变量 -> 逐层邻居哈希 -> 置位。返回 Python int 位集（bit_count 是 C 速度）。"""
    init=[]
    for k in range(len(atoms)):
        deg=len(adj.get(k,{}))
        init.append(hash((atoms[k][0], deg, atoms[k][1], sum(1 for b in adj.get(k,{}).values() if b==2))) & 0xffffffff)
    cur=init[:]
    bits=0
    for k in range(len(atoms)): bits |= 1 << (cur[k] % nbits)
    for _ in range(radius):
        nxt=[]
        for k in range(len(atoms)):
            nbrs=sorted(cur[j] for j in adj.get(k,{}))
            nxt.append(hash((cur[k], tuple(nbrs))) & 0xffffffff)
        cur=nxt
        for k in range(len(atoms)): bits |= 1 << (cur[k] % nbits)
    return bits

def tanimoto(a,b):
    u=(a|b).bit_count()
    return (a&b).bit_count()/u if u else 0.0
