import itertools, random, numpy as np, sympy as sp
from ref44 import *
from octalg import *
def octo_C(gamma=+1):
    tO={A:tau[A] for A in range(4)}
    for k in range(4,8): tO[k]=tau[k]*sigma
    if gamma<0: tO={A:tau[A] for A in range(8)}  # split (use 4+4 L with metric factor)
    return tO
# octonion structure constants
tO={A:tau[A] for A in range(4)}
for k in range(4,8): tO[k]=tau[k]*sigma
u=unit
E=sp.Matrix([list(u.T*tO[B].T) for B in range(8)]); Ei=E.T
LO={B: E*tO[B]*Ei for B in range(8)}
CO=mult_from_left({B: np.array(LO[B].tolist(),dtype=object) for B in range(8)})
LS=split_left_matrices()
CS=mult_from_left({B-1: np.array(LS[B].tolist(),dtype=object) for B in range(1,9)})
def qm(C,a,b):  # product inside Q=span(e0..e3) using table C
    x=list(a)+[0]*4; y=list(b)+[0]*4
    r=mul(C,x,y); assert r[4:]==[0]*4; return tuple(r[:4])
def qc(a): return (a[0],-a[1],-a[2],-a[3])
def qa(a,b): return tuple(x+y for x,y in zip(a,b))
def qs(a,b): return tuple(x-y for x,y in zip(a,b))
def ell_times(C,b):  # coefficients: ell*b = b0 e4 + b1 e5 + ... ? check
    return mul(C,[0,0,0,0,1,0,0,0],list(b)+[0]*4)
print('ell*e1..e3 (octonion):',[ell_times(CO,(0,1,0,0)),ell_times(CO,(0,0,1,0)),ell_times(CO,(0,0,0,1))])
print('e1..e3*ell (octonion):',[mul(CO,[0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0])])
cands={
 'Baez  (ac - d bbar, abar d + c b)': lambda C,a,b,c,d,g: (qs(qm(C,a,c),qm(C,d,qc(b))) if g<0 else qa(qm(C,a,c),qm(C,d,qc(b))), qa(qm(C,qc(a),d),qm(C,c,b))),
 'Wiki  (ac - dbar b, d a + b cbar)': lambda C,a,b,c,d,g: (qs(qm(C,a,c),qm(C,qc(d),b)) if g<0 else qa(qm(C,a,c),qm(C,qc(d),b)), qa(qm(C,d,a),qm(C,b,qc(c)))),
 'Wiki-opp (ac - b dbar, a d + cbar b)': lambda C,a,b,c,d,g: (qs(qm(C,a,c),qm(C,b,qc(d))) if g<0 else qa(qm(C,a,c),qm(C,b,qc(d))), qa(qm(C,a,d),qm(C,qc(c),b))),
 'Baez-opp (ac - bbar d, d abar + b c)': lambda C,a,b,c,d,g: (qs(qm(C,a,c),qm(C,qc(b),d)) if g<0 else qa(qm(C,a,c),qm(C,qc(b),d)), qa(qm(C,d,qc(a)),qm(C,b,c))),
}
rnd=random.Random(3)
R=lambda: tuple(rnd.randint(-4,4) for _ in range(4))
for idname,ident in [('x=a+ell*b',lambda a,b: list(a)+list(b)), ]:
    for tabname,C,g in [('octonion',CO,-1),('split',CS,+1)]:
        for cname,f in cands.items():
            ok=True
            for _ in range(200):
                a,b,c,d=R(),R(),R(),R()
                x=ident(a,b); y=ident(c,d)
                p,q=f(C,a,b,c,d,g)
                if mul(C,x,y)!=list(p)+list(q): ok=False;break
            print(f'{tabname:9s} {idname:10s} {cname:40s} match:{ok}')
