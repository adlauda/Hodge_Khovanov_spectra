"""Exact signed-cube, all-1 Tait graph, tree activity and orientation-cut check."""
from collections import Counter
from itertools import combinations
import json
from independent_kh import Khovanov, SparseMatrix, braid_closure


def orientation(s,word,reverse):
    labels=list(range(s)); events=[]
    for g in word:
        p=abs(g)-1
        events.append((labels[p],labels[p+1],1 if g>0 else -1))
        labels[p],labels[p+1]=labels[p+1],labels[p]
    comps={}
    for t in range(s):
        if t in comps: continue
        k=max(comps.values(),default=-1)+1
        while t not in comps:
            comps[t]=k
            t=labels[t]
    ell=max(comps.values())+1
    sums=Counter(); signs=[]
    for a,b,c in events:
        a,b=comps[a],comps[b]
        sign=c*(-1)**((a in reverse)+(b in reverse))
        signs.append(sign)
        if a!=b: sums[tuple(sorted([a,b]))]+=sign
    assert all(v%2==0 for v in sums.values())
    e=Counter()
    for S in range(1<<(ell-1)):
        mask=S<<1
        e[sum(v for (a,b),v in sums.items() if (mask>>a^mask>>b)&1)]+=1
    return ell,signs,sums,e


def table(kh, mark, reduced):
    mark=kh.marks[str(mark)]
    dims={}; ranks={}; m=kh.c
    def valid(state,labels):
        return bool((labels>>kh.resolution(state).vertex_circle[mark])&1)
    for i in range(-kh.diagram.nneg,kh.diagram.npos+1):
        for j in range(-3*m-4,3*m+5):
            source=kh.chain(i,j); target=kh.chain(i+1,j)
            if not len(source): continue
            d=kh.differential(i,j)
            if not reduced:
                dims[i,j]=len(source); ranks[i,j]=d.rank(); continue
            idx=[k for k,v in enumerate(source.basis) if valid(*v)]
            nxt={k for k,v in enumerate(target.basis) if valid(*v)}
            col=[d.columns[k] for k in idx]
            assert all(set(w).issubset(nxt) for w in col)
            dims[i,j]=len(idx); ranks[i,j]=SparseMatrix(len(target),col).rank()
    counts={p:v-ranks[p]-ranks.get((p[0]-1,p[1]),0) for p,v in dims.items()}
    assert all(v>=0 for v in counts.values())
    return Counter({p:v for p,v in counts.items() if v})


def spanning(kh):
    full=kh.resolution((1<<kh.c)-1)
    vs=full.count; m=kh.c
    edges=[(full.vertex_circle[cr[0]],full.vertex_circle[cr[2]]) for cr in kh.crossings]
    def tree(S):
        if len(S)!=vs-1: return False
        seen={0}
        for _ in range(vs):
            for k in S:
                a,b=edges[k]
                if a in seen or b in seen: seen.update((a,b))
        return len(seen)==vs
    ts={frozenset(T) for T in combinations(range(m),vs-1) if tree(T)}
    r=Counter()
    for T in ts:
        u=0
        for k in range(m):
            internal=k in T
            changers=[l for l in range(k) if (l in T)!=internal and T^{k,l} in ts]
            if not changers: u+=1 if internal else -1
        v=vs-1; w=kh.diagram.npos-kh.diagram.nneg; kap=m+2*v
        i2=2*u-4*v+w+kap; j2=4*u-4*v+3*w+kap-2
        assert i2%2==0 and j2%2==0
        r[i2//2,j2//2]+=1
    return r


def check(name,s,braid,reverse=()):
    ell,signs,lk2,e=orientation(s,braid,reverse)
    D=braid_closure(s,braid); D.npos=signs.count(1); D.nneg=signs.count(-1)
    kh=Khovanov(D)
    red=table(kh,0,True); ordinary=table(kh,0,False)
    st=spanning(kh)
    assert red==st,(name,red,st)
    diagonals={j-2*i for i,j in red}; assert len(diagonals)==1
    delta=diagonals.pop(); r=Counter({i:v for (i,j),v in st.items()})
    pred=Counter(); k=0
    for i in range(-kh.diagram.nneg-1,kh.diagram.npos+2):
        old=k; k=r[i]-e[i]-old; assert k>=0
        if e[i]+k: pred[i,2*i+delta]=e[i]+k
        if e[i]+old: pred[i,2*i+delta+2]=e[i]+old
    assert k==0 and ordinary==pred,(name,pred,ordinary)
    # Evaluate the two signed one-draw score numerators exactly over every
    # tree and orientation outcome. This is separate from the recurrence.
    for i in range(-kh.diagram.nneg-2,kh.diagram.npos+3):
        def score(offset,h):
            return (-1)**(i-offset-h) if h<=i-offset else 0
        low=sum(v*score(0,h) for h,v in r.items())+sum(v*score(1,h) for h,v in e.items())
        high=sum(v*score(1,h) for h,v in r.items())+sum(v*score(0,h) for h,v in e.items())
        assert (low,high)==(ordinary[i,2*i+delta],ordinary[i,2*i+delta+2])
    assert sum(ordinary.values())==sum(red.values())+2**(ell-1)
    def fmt(t): return [[i,j,k] for (i,j),k in sorted(t.items())]
    return dict(name=name,braid=braid,reverse=list(reverse),components=ell,
        linking={str(k):val//2 for k,val in sorted(lk2.items())},orientation_histogram=dict(sorted(e.items())),
        reduced=fmt(red),ordinary=fmt(ordinary),tree_count=sum(red.values()),ordinary_total=sum(ordinary.values()))


if __name__=="__main__":
    cases=[("unknot",1,[]),("positive curl",2,[1]),("negative curl",2,[-1]),
           ("positive Hopf",2,[1,1]),("positive Hopf, one component reversed",2,[1,1],(1,)),
           ("right trefoil",2,[1]*3),("left trefoil",2,[-1]*3),("figure eight",3,[1,-2]*2),
           ("Borromean rings",3,[1,-2]*3),("three-component Hopf chain",3,[1,1,-2,-2]),
           ("Hopf chain, one end reversed",3,[1,1,-2,-2],(0,))]
    data=[check(*item) for item in cases]
    print(json.dumps(data,indent=2))
