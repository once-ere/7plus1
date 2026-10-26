import numpy as np, sympy as sp, itertools
from ref71 import *
np.set_printoptions(linewidth=200)
T=[npm(T16O[A]) for A in range(8)]
E9=npm(T16O[8])
I16n=np.eye(16)
def realify(M):  # complex n x n -> real 2n x 2n
    X=M.real; Y=M.imag
    return np.block([[X,-Y],[Y,X]])
def commutant_basis(Gs):
    n=Gs[0].shape[0]
    rows=[]
    for G in Gs:
        # vec(MG - GM) = (G^T kron I - I kron G) vec(M)   (column-major vec)
        rows.append(np.kron(G.T,np.eye(n))-np.kron(np.eye(n),G))
    A=np.vstack(rows)
    u,s,vh=np.linalg.svd(A)
    tol=1e-9
    null=vh[np.sum(s>tol):]
    return [v.reshape(n,n,order='F') for v in null]
def classify(pos):   # pos: list of 8 booleans (True -> generator squares to +1)
    Gc=[ (T[A] if pos[A] else 1j*T[A]) for A in range(8)]
    # check
    for A in range(8):
        for B in range(8):
            target=2*(1 if pos[A] else -1)*(A==B)*np.eye(16)
            assert np.allclose(Gc[A]@Gc[B]+Gc[B]@Gc[A],target)
    Gr=[realify(G) for G in Gc]
    C=commutant_basis(Gr)
    dim=len(C)
    # orthonormalize w.r.t. trace form, remove identity component
    Id=np.eye(32)
    traceless=[]
    for M in C:
        M=M-np.trace(M)/32*Id
        traceless.append(M)
    # basis of traceless part (rank 3)
    V=np.array([m.flatten() for m in traceless])
    u,s,vh=np.linalg.svd(V)
    r=np.sum(s>1e-9)
    Bs=[vh[i].reshape(32,32) for i in range(r)]
    Gram=np.array([[np.trace(a@b)/32 for b in Bs] for a in Bs])
    ev=np.linalg.eigvalsh((Gram+Gram.T)/2)
    return dim,r,np.round(ev,6)
for p in range(8,-1,-1):
    q=8-p
    pos=[True]*p+[False]*q
    dim,r,ev=classify(pos)
    kind='H-type (division) => H(8), no real 16-dim rep' if all(ev<0) else 'R-type => R(16)'
    print(f'(p,q)=({p},{q}) commutant dim {dim}, traceless rank {r}, Gram eigenvalues {ev}: {kind}')
