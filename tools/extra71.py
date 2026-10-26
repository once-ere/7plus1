"""Extra independent checks (Python) of the claims used in 7+1.nb."""
import itertools, numpy as np, sympy as sp
from ref71 import *
R={}
# vielbein relation
R['EO == eta4488 . E44[8]']=sp.simplify(EAaO-eta4488*EAa(8))==sp.zeros(8)
# Gamma53 from 4+4 (split trace form)
G53={0:T16A[0]}
for i in range(1,8): G53[i]=T16A[0]*T16A[i]
eta53=sp.diag(1,-1,-1,-1,1,1,1,1)
R['G53 Cl(5,3) real']=anticomm_ok(G53,eta53,16)
# Gamma08
G08={A:T16O[8]*T16O[A] for A in range(8)}
R['G08 Cl(0,8) real']=anticomm_ok(G08,-sp.eye(8),16)
# so(1,7) relations
G=G17; eta=eta17
Sg=lambda A,B:(G[A]*G[B]-G[B]*G[A])/4
okS=True
for A,B,C,Dd in itertools.product(range(8),repeat=4):
    if A>=B or C>=Dd: continue
    lhs=Sg(A,B)*Sg(C,Dd)-Sg(C,Dd)*Sg(A,B)
    rhs=eta[B,C]*Sg(A,Dd)-eta[A,C]*Sg(B,Dd)-eta[B,Dd]*Sg(A,C)+eta[A,Dd]*Sg(B,C)
    if lhs!=rhs: okS=False;break
R['so(1,7): [S_AB,S_CD]=eta_BC S_AD - eta_AC S_BD - eta_BD S_AC + eta_AD S_BC']=okS
okG=all(Sg(A,B)*G[C]-G[C]*Sg(A,B)==eta[B,C]*G[A]-eta[A,C]*G[B] for A in range(8) for B in range(8) for C in range(8))
R['[S_AB,G_C]=eta_BC G_A - eta_AC G_B']=okG
# boost
b=sp.symbols('beta',real=True)
S=sp.cosh(b/2)*I16+sp.sinh(b/2)*G[0]*G[1]
Sinv=sp.cosh(b/2)*I16-sp.sinh(b/2)*G[0]*G[1]
R['S Sinv = 1']=sp.simplify(S*Sinv-I16)==sp.zeros(16)
# Lambda from traces: S G_C S^-1 = sum_B L[B,C] G_B ; L[B,C] = Tr(G^B S G_C S^-1)/16 with G^B = eta^BB G_B
L=sp.zeros(8)
for B_ in range(8):
    for C_ in range(8):
        L[B_,C_]=sp.simplify((eta[B_,B_]*G[B_]*S*G[C_]*Sinv).trace()/16)
print('Lambda boost (rows B = 0..7):')
for r_ in range(8): print('   ', [str(L[r_, c_]) for c_ in range(8)])
R['S G_C S^-1 == sum_B L[B,C] G_B']=all(sp.simplify(S*G[C_]*Sinv-sum((L[B_,C_]*G[B_] for B_ in range(8)),sp.zeros(16)))==sp.zeros(16) for C_ in range(8))
R['L^T eta L == eta']=sp.simplify(L.T*eta*L-eta)==sp.zeros(8)
# rotation 4-5
th=sp.symbols('theta',real=True)
Rr=sp.cos(th/2)*I16+sp.sin(th/2)*G[4]*G[5]
Rinv=sp.cos(th/2)*I16-sp.sin(th/2)*G[4]*G[5]
Lr=sp.zeros(8)
for B_ in range(8):
    for C_ in range(8):
        Lr[B_,C_]=sp.simplify((eta[B_,B_]*G[B_]*Rr*G[C_]*Rinv).trace()/16)
print('Lambda rotation (rows B = 0..7):')
for r_ in range(8): print('   ', [str(Lr[r_, c_]) for c_ in range(8)])
# on-shell example
p=[3,2,2,0,0,0,0,0]
pv=sp.Matrix(p)
R['p on shell p.eta17.p == 1']=(pv.T*eta17*pv)[0]==1
ps=sum((p[A]*G[A] for A in range(8)),sp.zeros(16))
R['pslash^2 == 1']=ps*ps==I16
ns=(ps-I16).nullspace()
R['8 real solutions']=len(ns)==8 and all(all(x.is_rational for x in v) for v in ns)
# Frobenius-Schur type indicator for all signatures
def intertwiner(Gs):
    n=16
    rows=[np.kron(np.eye(n),np.conj(Gm).T)-np.kron(Gm,np.eye(n)) for Gm in Gs]   # row-major vec: vec(M Gbar)= (1 kron Gbar^T) vec M ; vec(G M) = (G kron 1) vec M
    A=np.vstack(rows)
    u,s,vh=np.linalg.svd(A)
    null=vh[np.sum(s>1e-9):]
    return [np.conj(v).reshape(n,n) for v in null]
T=[np.array(T16O[A].tolist(),dtype=complex) for A in range(8)]
types=[]
for pp in range(8,-1,-1):
    Gs=[T[A] if A<pp else 1j*T[A] for A in range(8)]
    for A in range(8):
        for B_ in range(8):
            assert np.allclose(Gs[A]@Gs[B_]+Gs[B_]@Gs[A],2*(1 if A<pp else -1)*(A==B_)*np.eye(16))
    Ms=intertwiner(Gs)
    M=Ms[0]
    # verify intertwining
    assert all(np.allclose(M@np.conj(Gm),Gm@M) for Gm in Gs)
    MM=M@np.conj(M)
    lam=MM[0,0]
    assert np.allclose(MM,lam*np.eye(16))
    types.append(((pp,8-pp),len(Ms),float(round(lam.real,6)),'R(16)' if lam.real>0 else 'H(8)'))
print('indicator classification (signature, number of intertwiners, M.conj(M) for the normalized M, algebra):')
for t_ in types: print('    Cl%s: %d intertwiner, M.conj(M) = %+.4f * 1  ->  %s' % (str(t_[0]).replace(' ', ''), t_[1], t_[2], t_[3]))
R['classification']=[t[3] for t in types]==['R(16)','H(8)','H(8)','R(16)','R(16)','H(8)','H(8)','R(16)','R(16)']
# commutant of real 32x32 G71R32: explicit quaternion units
s1=sp.Matrix([[0,1],[1,0]]); s3=sp.Matrix([[1,0],[0,-1]]); e2=sp.Matrix([[0,1],[-1,0]])
G8=G17[8]
Qi=sp.kronecker_product(e2,I16); Qj=sp.kronecker_product(s1,G8); Qk=sp.kronecker_product(s3,G8)
comm=all(Q*G71R32[A]==G71R32[A]*Q for Q in (Qi,Qj,Qk) for A in range(8))
quat=(Qi*Qi==-sp.eye(32)) and (Qj*Qj==-sp.eye(32)) and (Qk*Qk==-sp.eye(32)) and (Qi*Qj==Qk) and (Qj*Qi==-Qk)
R['explicit quaternion commutant of G71R32']=comm and quat
# octonion basis form of G17
EO16=blk(EAaO,Z8,Z8,EAaO)
G17e={A: sp.simplify(EO16*G17[A]*EO16.T) for A in range(8)}
R['G17e[0]==[[0,1],[1,0]]x1, G17e[i]==diag(LO[i],-LO[i])']=G17e[0]==blk(Z8,I8,I8,Z8) and all(G17e[i]==blk(LO[i],Z8,Z8,-LO[i]) for i in range(1,8))
xs=sp.symbols('x0:8')
Gx=sum((xs[A]*G17[A] for A in range(8)),sp.zeros(16))
R['Gamma17(x)^2 == (x.eta17.x) 1']=sp.expand(Gx*Gx-(sum(eta17[A,A]*xs[A]**2 for A in range(8)))*I16)==sp.zeros(16)
# ---- real spinor fields (section 5.5 of 7+1.nb) ----
# (a) the matrices anticommuting with every G17[A]: numerically a 1-dimensional space, containing G17[8]
Gn=[np.array(G17[A].tolist(),dtype=float) for A in range(8)]
sysA=np.vstack([np.kron(np.eye(16),Gm.T)+np.kron(Gm,np.eye(16)) for Gm in Gn])
R['anticommutant of G17: dimension 1 (numerical rank 255 of 2048x256 system)']=int(np.linalg.matrix_rank(sysA))==255
R['G17[8] anticommutes with all G17[A], G17[8]^2 == -1']=all(G17[8]*G17[A]==-G17[A]*G17[8] for A in range(8)) and G17[8]*G17[8]==-I16
# (b) a REAL massless field: k = (1,1,0,...,0), Psi = w cos(k.x)
X_=sp.Matrix(xs)
def dirac(G,metric,Psi):
    return sum((metric[A,A]*G[A]*Psi.diff(xs[A]) for A in range(8)),sp.zeros(Psi.shape[0],1))
kv=sp.Matrix([1,1,0,0,0,0,0,0])
ks=sum((kv[A]*G17[A] for A in range(8)),sp.zeros(16))
nk=ks.nullspace()
R['null k: k.eta17.k == 0, kslash^2 == 0, 8 rational solutions of kslash w = 0']=(kv.T*eta17*kv)[0]==0 and ks*ks==sp.zeros(16) and len(nk)==8 and all(all(x.is_rational for x in v) for v in nk)
eta17inv=eta17.inv()
Psi0=nk[0]*sp.cos((kv.T*eta17*X_)[0])
R['Psi = w cos(k.x) is a REAL solution of G17^A d_A Psi = 0']=sp.simplify(dirac(G17,eta17inv,Psi0))==sp.zeros(16,1)
# (c) a REAL massive field with the real 32x32 generators: Psi = w cos(p.x) + (P.w/m) sin(p.x)
eta71inv=eta71.inv()
P32=sum((p[A]*G71R32[A] for A in range(8)),sp.zeros(32))
R['p.eta71.p == -m^2 and P32^2 == -m^2']=(pv.T*eta71*pv)[0]==-1 and P32*P32==-sp.eye(32)
wv=sp.Matrix(sp.symbols('w1:33',real=True))
ph=(pv.T*eta17*X_)[0]
Psi32=wv*sp.cos(ph)+P32*wv*sp.sin(ph)
R['Psi32 (32 arbitrary real w) solves (G71R32^A d_A - m) Psi = 0']=sp.simplify(dirac(G71R32,eta71inv,Psi32)-Psi32)==sp.zeros(32,1)
wave17=lambda Psi: sum((eta17inv[A,A]*Psi.diff(xs[A],2) for A in range(8)),sp.zeros(Psi.shape[0],1))
R['Psi32 obeys the Klein-Gordon equation (box17 + m^2) Psi = 0']=sp.simplify(wave17(Psi32)+Psi32)==sp.zeros(32,1)
fs=sp.Matrix([sp.Function('f%d'%k)(*xs) for k in range(32)])
R['(G71R32^A d_A)^2 == -box17 on 32 arbitrary functions']=sp.expand(dirac(G71R32,eta71inv,dirac(G71R32,eta71inv,fs))+wave17(fs))==sp.zeros(32,1)
for k,v in R.items(): print(f'{k:70s}: {v}')
