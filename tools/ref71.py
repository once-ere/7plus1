"""Reference (Python/sympy/numpy) implementation of everything planned for 7+1.nb.
All integer matrices are kept exact (sympy) where it matters; heavy linear algebra uses numpy."""
import itertools, random
import numpy as np
import sympy as sp
from ref44 import *            # 4+4.nb objects: tau, taubar, sigma, T16A, eta4488, unit, U, ...
from octalg import *

I8=sp.eye(8); I16=sp.eye(16)
def npm(M, dtype=float): return np.array(M.tolist(), dtype=dtype)

# ---------------- octonion 8x8 matrices (mirror of tau) ----------------
tauO={A: tau[A] for A in range(4)}
for k in range(4,8): tauO[k]=tau[k]*sigma
tauObar={A: tauO[A].T for A in range(8)}
gO=sp.eye(8)                     # octonion norm metric (mirror of g_AB = eta4488)

# vielbein with the SAME unit spinor as 4+4.nb (hUSE = 8)
uO=unit
EAaO=sp.Matrix([list(uO.T*tauObar[B]) for B in range(8)])     # rows B
EaAO=sp.simplify(EAaO.inv())
assert EaAO==sp.Matrix.hstack(*[tauO[B]*uO for B in range(8)])
LO={B: sp.simplify(gO[B,B]*EAaO*tauO[B]*EaAO) for B in range(8)}   # left-mult matrices
for B in range(8): assert all(x in (0,1,-1) for x in LO[B])
CO=mult_from_left({B: np.array(LO[B].tolist(),dtype=object) for B in range(8)})
# split (4+4) structure constants, 0-based relabel eA[B] -> eS[B-1]
LS44=split_left_matrices()
LS={B-1: LS44[B] for B in range(1,9)}
CS=mult_from_left({B: np.array(LS[B].tolist(),dtype=object) for B in range(8)})

def omul(x,y): return mul(CO,x,y)
def smul(x,y): return mul(CS,x,y)

# ---------------- Clifford: Euclidean mirror T16O ----------------
def blk(a,b,c,d): return sp.Matrix(sp.BlockMatrix([[a,b],[c,d]]))
Z8=sp.zeros(8)
T16O={A: blk(Z8,tauObar[A],tauO[A],Z8) for A in range(8)}
T16O[8]=I16
for A in range(8): T16O[8]=T16O[8]*T16O[A]
S44=T16A[4]*T16A[5]*T16A[6]*T16A[7]
# ---------------- 7+1 real (Cl_{1,7}) ----------------
# eta17 = diag(+1,-1,...,-1) is simply the squares of the octonion basis: e_0^2 = +1 (the time direction x0),
# e_i^2 = -1 (the seven space directions), eta17_AB = Re(e_A e_B).  The spacetime has one time-like and seven
# space-like dimensions; the tensor package xAct (section 12 of 7+1.nb) declares it with the signature list {7, 1, 0}:
# seven space-like, one time-like, none degenerate -- the counts of +1 and -1 in eta71 = -eta17 (4+4.nb convention).
eta17=sp.diag(1,-1,-1,-1,-1,-1,-1,-1)
eta71=-eta17
G17={0: T16O[0]}
for i in range(1,8): G17[i]=T16O[0]*T16O[i]
G17[8]=I16
for A in range(8): G17[8]=G17[8]*G17[A]
# ---------------- 7+1 complex (Cl_{7,1}, 4+4 sign convention) ----------------
tauC={0: sp.I*I8}; tauCbar={0: sp.I*I8}
for k in range(1,8): tauC[k]=tauO[k]; tauCbar[k]=tauObar[k]
G71={A: blk(Z8,tauCbar[A],tauC[A],Z8) for A in range(8)}
G71[8]=I16
for A in range(8): G71[8]=G71[8]*G71[A]
Ssim=(I16-sp.I*G17[0])/sp.sqrt(2)
# ---------------- 7+1 real 32x32 ----------------
eps2=sp.Matrix([[0,1],[-1,0]])
G71R32={A: sp.kronecker_product(eps2,G17[A]) for A in range(8)}

def anticomm_ok(G,metric,n):
    Id=sp.eye(n)
    return all(sp.simplify(G[A]*G[B]+G[B]*G[A]-2*metric[A,B]*Id)==sp.zeros(n) for A in range(8) for B in range(8))

def products(G):
    """all 256 ordered products G[a1]...G[ak], a1<...<ak (k=0..8)"""
    out=[]
    n=G[0].shape[0]
    for k in range(9):
        for S in itertools.combinations(range(8),k):
            M=sp.eye(n)
            for a in S: M=M*G[a]
            out.append((S,M))
    return out

if __name__=='__main__':
    R={}
    # octonion matrices
    R['tauO antisymmetric (A>=1)']=all(tauO[A].T==-tauO[A] for A in range(1,8))
    R['tauO[A]^2=-1 (A>=1)']=all(tauO[A]*tauO[A]==-I8 for A in range(1,8))
    R['tauO Euclid relation']=all(tauO[A]*tauObar[B]+tauO[B]*tauObar[A]==2*gO[A,B]*I8 for A in range(8) for B in range(8))
    R['tauObar tauO relation']=all(tauObar[A]*tauO[B]+tauObar[B]*tauO[A]==2*gO[A,B]*I8 for A in range(8) for B in range(8))
    P=I8
    for A in range(1,8): P=P*tauO[A]
    R['tauO1..7 product']=list(P.diagonal())
    R['sigma commutes tau1..3, anticommutes tau4..7']=all(sigma*tau[h]==tau[h]*sigma for h in range(4)) and all(sigma*tau[k]==-tau[k]*sigma for k in range(4,8))
    R['E orthogonal']=EAaO*EAaO.T==I8
    res=report(CO,[1]*8,'oct')
    R.update({'oct: '+k:v for k,v in res.items()})
    # sign relation octonion vs split
    s=lambda A,B: -1 if (A>=4 and B>=4) else 1
    R['e_A*e_B (O) = s_AB e_A*e_B (S)']=all(omul(basis(A),basis(B))==scal(s(A,B),smul(basis(A),basis(B))) for A in range(8) for B in range(8))
    # Clifford
    R['T16O Euclid Cl(8,0)']=anticomm_ok(T16O,sp.eye(8),16)
    R['T16O[A]==T16A44[A] (A<=3)']=all(T16O[A]==T16A[A] for A in range(4))
    R['T16O[k]==T16A44[k].S44 (k>=4)']=all(T16O[k]==T16A[k]*S44 for k in range(4,8))
    R['S44==diag(sigma,sigma)']=S44==blk(sigma,Z8,Z8,sigma)
    R['T16O[8] diag']=list(T16O[8].diagonal())
    R['eta17 = squares of the octonion basis: Re(e_A e_B) = eta17_AB']=sp.Matrix(8,8,lambda A,B: omul(basis(A),basis(B))[0])==eta17
    R['G17 Cl(1,7)']=anticomm_ok(G17,eta17,16)
    R['G17 real']=all(all(x.is_integer for x in G17[A]) for A in range(8))
    R['G17[0] sym, G17[i] antisym']=G17[0].T==G17[0] and all(G17[i].T==-G17[i] for i in range(1,8))
    R['G17[8]^2']=(G17[8]*G17[8])[0,0]
    R['G17 blocks: G17[i]==diag(tauO,tauObar)']=all(G17[i]==blk(tauO[i],Z8,Z8,tauObar[i]) for i in range(1,8))
    R['G71 Cl(7,1) eta71']=anticomm_ok(G71,eta71,16)
    R['G71[8]^2']=sp.simplify((G71[8]*G71[8])[0,0])
    R['G71 = Ssim (I G17) Ssim^-1']=all(sp.simplify(Ssim*(sp.I*G17[A])*Ssim.inv()-G71[A])==sp.zeros(16) for A in range(8))
    R['G71R32 Cl(7,1) real']=anticomm_ok(G71R32,eta71,32)
    for k,v in R.items(): print(f'{k:50s}: {v}')
