"""Exact reproduction of the 4+4.nb definitions (sympy, exact arithmetic)."""
import itertools
import sympy as sp

def signature(perm):
    """Mathematica Signature[] of a list (0 if repeated)."""
    perm=list(perm)
    if len(set(perm))<len(perm): return 0
    s=1
    p=perm[:]
    for i in range(len(p)):
        for j in range(i+1,len(p)):
            if p[i]>p[j]: s=-s
    return s

def kd(a,b): return 1 if a==b else 0

ID4=sp.eye(4); ID8=sp.eye(8); ID16=sp.eye(16)
eta4488=sp.diag(1,1,1,1,-1,-1,-1,-1)
Z4=sp.zeros(4); Z8=sp.zeros(8)

def blocks(a,b,c,d):
    return sp.Matrix(sp.BlockMatrix([[a,b],[c,d]]))

sigma=blocks(Z4,ID4,ID4,Z4)          # σ = σ88
G16=sp.diag(eta4488,-eta4488)

def Qa(h,p,q): return signature([h,p,q,4])
def Qb(h,p,q): return ID4[p-1,3]*ID4[q-1,h-1]-ID4[p-1,h-1]*ID4[q-1,3]
def SD(h,p,q): return Qa(h,p,q)-Qb(h,p,q)
def ASD(h,p,q): return Qa(h,p,q)+Qb(h,p,q)

s4={h: sp.Matrix(4,4,lambda p,q: SD(h,p+1,q+1)) for h in (1,2,3)}
t4={h: sp.Matrix(4,4,lambda p,q: ASD(h,p+1,q+1)) for h in (1,2,3)}

zeta={0:ID8}
for h in (1,2,3):
    zeta[7-h]=blocks(Z4,t4[h],-t4[h],Z4)
for h in (1,2,3):
    zeta[h]=blocks(Z4,s4[h],s4[h],Z4)
zeta[7]=zeta[1]*zeta[2]*zeta[3]*zeta[4]*zeta[5]*zeta[6]

tau={A:zeta[A] for A in range(8)}
taubar={0:ID8}
for A in range(1,8):
    taubar[A]=sigma*tau[A].T*sigma

def T16(A):
    return blocks(Z8,taubar[A],tau[A],Z8)
T16A={A:T16(A) for A in range(8)}
T16A[8]=sp.eye(16)
for A in range(8): T16A[8]=T16A[8]*T16A[A]
sigma16=T16A[0]*T16A[1]*T16A[2]*T16A[3]

# eigen-basis of σ as in 4+4.nb (Eigensystem order):
evecs=[ [0,0,0,-1,0,0,0,1],[0,0,-1,0,0,0,1,0],[0,-1,0,0,0,1,0,0],[-1,0,0,0,1,0,0,0],
        [0,0,0,1,0,0,0,1],[0,0,1,0,0,0,1,0],[0,1,0,0,0,1,0,0],[1,0,0,0,1,0,0,0]]
U=[sp.Matrix(v)/sp.sqrt(2) for v in evecs]   # 𝕌[[h]] , h=1..8 -> U[h-1]
hUSE=8
unit=U[hUSE-1]

def EAa(h):
    """Subscript[(E)^(A), a][h] = Table[𝕌[[h]].σ.τ̄[B], {B,0,7}]  (8x8, row B)"""
    u=U[h-1]
    return sp.Matrix([ list((u.T*sigma*taubar[B])) for B in range(8)])
def EaA(h):
    return sp.simplify(EAa(h).inv())

gAB=eta4488

def split_left_matrices(h=hUSE):
    """L[B] (B=1..8, i.e. eA[B]) = g_AB[[B,B]] E^A_a[h].τ[B-1].E^a_A[h]"""
    E=EAa(h); Ei=EaA(h)
    return {B: sp.simplify(gAB[B-1,B-1]*E*tau[B-1]*Ei) for B in range(1,9)}

if __name__=='__main__':
    # checks from 4+4.nb
    print('ζ7 =', list(zeta[7].diagonal()))
    print('σ == ζ1ζ2ζ3 :', sigma==zeta[1]*zeta[2]*zeta[3])
    ok=all(((zeta[A]*zeta[B]+zeta[B]*zeta[A])/2)==-eta4488[A,B]*ID8 for A in range(1,8) for B in range(1,8))
    print('{ζA,ζB}/2 == -η ID8 (A,B=1..7):',ok)
    print('-η_AA ζA == ζA^T :',all(-eta4488[A,A]*zeta[A]==zeta[A].T for A in range(1,8)))
    ok=all(((tau[A]*taubar[B]+tau[B]*taubar[A])/2)==eta4488[A,B]*ID8 for A in range(8) for B in range(8))
    print('(τA τ̄B+τB τ̄A)/2 == η ID8 :',ok)
    ok=all(((taubar[A]*tau[B]+taubar[B]*tau[A])/2)==eta4488[A,B]*ID8 for A in range(8) for B in range(8))
    print('(τ̄A τB+τ̄B τA)/2 == η ID8 :',ok)
    print('τ1..τ7 product == ID8:', (tau[1]*tau[2]*tau[3]*tau[4]*tau[5]*tau[6]*tau[7])==ID8)
    print('σ == τ4τ5τ6τ7:', sigma==tau[4]*tau[5]*tau[6]*tau[7])
    ok=all(((T16A[A]*T16A[B]+T16A[B]*T16A[A])/2)==eta4488[A,B]*ID16 for A in range(8) for B in range(8))
    print('{T16A,T16B}/2 == η ID16 :',ok)
    print('T16[8] diag:', list(T16A[8].diagonal()))
    print('σ16 == diag(-σ,σ):', sigma16==blocks(-sigma,Z8,Z8,sigma))
    print('unit:', list(unit), ' unit.σ.unit =', sp.simplify((unit.T*sigma*unit)[0]))
    # E check as in notebook
    for h in range(1,9):
        u=U[h-1]
        lhs=EaA(h)
        rhs=sp.simplify((u.T*sigma*u)[0])*sp.Matrix.hstack(*[eta4488[B,B]*tau[B]*u for B in range(8)])
        assert sp.simplify(lhs-rhs)==sp.zeros(8), h
    print('E^a_A[h] == (u σ u) Transpose[Table[η_BB τ[B-1].u]] for all h: True')
    L=split_left_matrices()
    names=['eA[%d]'%k for k in range(1,9)]
    def rowtxt(B):
        out=[]
        for k in range(8):
            col=L[B][:,k]
            terms=[]
            for D in range(8):
                c=col[D]
                if c!=0: terms.append(('' if c==1 else '-' if c==-1 else str(c)+'*')+names[D])
            out.append('+'.join(terms))
        return out
    for B in range(1,9):
        print(names[B-1],':',rowtxt(B))
