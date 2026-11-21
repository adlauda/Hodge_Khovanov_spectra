"""Exact ordinary and one-basepoint reduced Kunneth check on split projections."""
from collections import Counter
from check_alternating import table
from independent_kh import Khovanov, braid_closure
import json

def tensor(tables):
    product=Counter({(0,0):1})
    for table in tables:
        updated=Counter()
        for (i,j),a in product.items():
            for (h,k),b in table.items():
                updated[i+h,j+k]+=a*b
        product=updated
    return product

def check(label,components):
    total_strands=0
    word=[]
    blocks=[]
    for s,w in components:
        word += [(1 if x>0 else -1)*(abs(x)+total_strands) for x in w]
        blocks.append(Khovanov(braid_closure(s,w)))
        total_strands+=s
    joined=Khovanov(braid_closure(total_strands,word))
    ordinary=table(joined,0,False); reduced=table(joined,0,True)
    assert ordinary==tensor([table(b,0,False) for b in blocks])
    assert reduced==tensor([table(b,0,k==0) for k,b in enumerate(blocks)])
    return dict(name=label,ordinary=sum(ordinary.values()),reduced=sum(reduced.values()),
                ordinary_bins=len(ordinary),reduced_bins=len(reduced))
if __name__=="__main__":
    U=(1,[]); H=(2,[1,1]); T=(2,[1,1,1])
    print(json.dumps([check(name,parts) for name,parts in
                      [("unlink", [U,U]),("Hopf plus unknot",[H,U]),("trefoil plus Hopf",[T,H])]],indent=2))
