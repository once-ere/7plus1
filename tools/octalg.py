"""Generic checks for an 8-dim algebra given by structure constants C[B][k] = vector (coeffs of e_D)
product: (x*y)_D = sum_{B,k} x_B y_k C[B][k][D]   (C[B][k] = coefficients of e_B e_k)."""
import itertools, random
import numpy as np

def mult_from_left(L):
    """L: dict B(0-based)-> 8x8 integer numpy array; e_B e_k = sum_D L[B][D,k] e_D"""
    C=np.zeros((8,8,8),dtype=object)
    for B in range(8):
        for k in range(8):
            for D in range(8):
                C[B,k,D]=int(L[B][D,k])
    return C

def mul(C,x,y):
    out=[0]*8
    for B in range(8):
        if x[B]==0: continue
        for k in range(8):
            if y[k]==0: continue
            for D in range(8):
                c=C[B,k,D]
                if c: out[D]+=x[B]*y[k]*c
    return out

def add(x,y): return [a+b for a,b in zip(x,y)]
def sub(x,y): return [a-b for a,b in zip(x,y)]
def scal(s,x): return [s*a for a in x]
def conj(x): return [x[0]]+[-a for a in x[1:]]

def basis(i):
    v=[0]*8; v[i]=1; return v

def assoc(C,x,y,z):
    return sub(mul(C,mul(C,x,y),z),mul(C,x,mul(C,y,z)))

def report(C,metric,name,ntest=300,seed=1):
    rnd=random.Random(seed)
    R=lambda: [rnd.randint(-5,5) for _ in range(8)]
    one=basis(0)
    res={}
    res['unit two-sided']=all(mul(C,one,basis(i))==basis(i) and mul(C,basis(i),one)==basis(i) for i in range(8))
    res['e_i^2']=[mul(C,basis(i),basis(i))[0] for i in range(1,8)]
    N=lambda x: sum(metric[i]*x[i]*x[i] for i in range(8))
    res['x xbar = N(x) 1']=all(mul(C,x,conj(x))==[N(x)]+[0]*7 and mul(C,conj(x),x)==[N(x)]+[0]*7 for x in [R() for _ in range(ntest)])
    ok=True
    for _ in range(ntest):
        x,y=R(),R()
        if N(mul(C,x,y))!=N(x)*N(y): ok=False;break
    res['N(xy)=N(x)N(y)']=ok
    ok=True
    for _ in range(ntest):
        x,y=R(),R()
        if conj(mul(C,x,y))!=mul(C,conj(y),conj(x)): ok=False;break
    res['conj(xy)=conj(y)conj(x)']=ok
    ok=True
    for _ in range(ntest):
        x,y=R(),R()
        if assoc(C,x,x,y)!=[0]*8 or assoc(C,y,x,x)!=[0]*8 or assoc(C,x,y,x)!=[0]*8: ok=False;break
    res['alternative (left,right,flexible)']=ok
    ok=True
    for _ in range(ntest//3):
        x,y,z=R(),R(),R()
        # Moufang: z(x(zy)) = ((zx)z)y ; x(z(yz)) = ((xz)y)z ; (zx)(yz) = (z(xy))z
        if mul(C,z,mul(C,x,mul(C,z,y)))!=mul(C,mul(C,mul(C,z,x),z),y): ok=False;break
        if mul(C,x,mul(C,z,mul(C,y,z)))!=mul(C,mul(C,mul(C,x,z),y),z): ok=False;break
        if mul(C,mul(C,z,x),mul(C,y,z))!=mul(C,mul(C,z,mul(C,x,y)),z): ok=False;break
    res['Moufang identities']=ok
    # non-associativity count on basis
    na=sum(1 for i,j,k in itertools.product(range(8),repeat=3) if assoc(C,basis(i),basis(j),basis(k))!=[0]*8)
    res['#nonassociative basis triples']=na
    nc=sum(1 for i,j in itertools.product(range(8),repeat=2) if mul(C,basis(i),basis(j))!=mul(C,basis(j),basis(i)))
    res['#noncommuting ordered basis pairs']=nc
    return res

def fano_triples(C):
    """triples (i,j,k) with e_i e_j = +e_k, i<j , i,j,k in 1..7"""
    T=[]
    for i in range(1,8):
        for j in range(1,8):
            if i==j: continue
            p=mul(C,basis(i),basis(j))
            k=[d for d in range(8) if p[d]!=0]
            assert len(k)==1
            k=k[0]
            if p[k]==1 and i<j: T.append((i,j,k))
    return T
