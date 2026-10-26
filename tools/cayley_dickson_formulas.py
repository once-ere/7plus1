import random
from fractions import Fraction
def qmul(p,q):
    a1,b1,c1,d1=p; a2,b2,c2,d2=q   # Hamilton: i j = k
    return (a1*a2-b1*b2-c1*c2-d1*d2,
            a1*b2+b1*a2+c1*d2-d1*c2,
            a1*c2-b1*d2+c1*a2+d1*b2,
            a1*d2+b1*c2-c1*b2+d1*a2)
def qconj(p): return (p[0],-p[1],-p[2],-p[3])
def qadd(p,q): return tuple(x+y for x,y in zip(p,q))
def qsub(p,q): return tuple(x-y for x,y in zip(p,q))
formulas={
 'user2: (ac - dbar b, abar d + c b)': lambda a,b,c,d: (qsub(qmul(a,c),qmul(qconj(d),b)), qadd(qmul(qconj(a),d),qmul(c,b))),
 'Baez: (ac - d bbar, abar d + c b)': lambda a,b,c,d: (qsub(qmul(a,c),qmul(d,qconj(b))), qadd(qmul(qconj(a),d),qmul(c,b))),
 'Wiki: (ac - dbar b, d a + b cbar)': lambda a,b,c,d: (qsub(qmul(a,c),qmul(qconj(d),b)), qadd(qmul(d,a),qmul(b,qconj(c)))),
 'user1-fixed: (ac - dbar b, d a + c bbar)': lambda a,b,c,d: (qsub(qmul(a,c),qmul(qconj(d),b)), qadd(qmul(d,a),qmul(c,qconj(b)))),
 'split(Wiki,gamma=+1): (ac + dbar b, d a + b cbar)': lambda a,b,c,d: (qadd(qmul(a,c),qmul(qconj(d),b)), qadd(qmul(d,a),qmul(b,qconj(c)))),
}
def omul(f,x,y):
    a,b=tuple(x[:4]),tuple(x[4:]); c,d=tuple(y[:4]),tuple(y[4:])
    r1,r2=f(a,b,c,d); return list(r1)+list(r2)
def sub(x,y): return [p-q for p,q in zip(x,y)]
rnd=random.Random(5)
R=lambda: [rnd.randint(-4,4) for _ in range(8)]
for name,f in formulas.items():
    alt=all(sub(omul(f,omul(f,x,x),y),omul(f,x,omul(f,x,y)))==[0]*8 and sub(omul(f,omul(f,y,x),x),omul(f,y,omul(f,x,x)))==[0]*8 for x,y in [(R(),R()) for _ in range(200)])
    N=lambda x: sum(t*t for t in x)
    comp=all(N(omul(f,x,y))==N(x)*N(y) for x,y in [(R(),R()) for _ in range(200)])
    Ns=lambda x: sum(t*t for t in x[:4])-sum(t*t for t in x[4:])
    comps=all(Ns(omul(f,x,y))==Ns(x)*Ns(y) for x,y in [(R(),R()) for _ in range(200)])
    unit=[1,0,0,0,0,0,0,0]
    u=all(omul(f,unit,[int(i==k) for i in range(8)])==[int(i==k) for i in range(8)] and omul(f,[int(i==k) for i in range(8)],unit)==[int(i==k) for i in range(8)] for k in range(8))
    print(f'{name:55s} unit:{u} alternative:{alt} N_euclid multiplicative:{comp} N_split multiplicative:{comps}')
print('--- with opposite quaternion product (ij=-k) ---')
_q=qmul
qmul=lambda p,q: _q(q,p)
for name,f in formulas.items():
    alt=all(sub(omul(f,omul(f,x,x),y),omul(f,x,omul(f,x,y)))==[0]*8 and sub(omul(f,omul(f,y,x),x),omul(f,y,omul(f,x,x)))==[0]*8 for x,y in [(R(),R()) for _ in range(200)])
    N=lambda x: sum(t*t for t in x)
    comp=all(N(omul(f,x,y))==N(x)*N(y) for x,y in [(R(),R()) for _ in range(200)])
    print(f'{name:55s} alternative:{alt} N_euclid multiplicative:{comp}')
# explicit counterexample for user2 with Hamilton product
qmul=_q
f=formulas['user2: (ac - dbar b, abar d + c b)']
import itertools
B=lambda i: [int(k==i) for k in range(8)]
found=None
for i,j in itertools.product(range(8),repeat=2):
    x,y=B(i),B(j)
    lhs=omul(f,omul(f,x,x),y); rhs=omul(f,x,omul(f,x,y))
    if lhs!=rhs: found=(i,j,lhs,rhs);break
print('basis counterexample to x^2 y = x(xy):',found)
for i,j in itertools.product(range(8),repeat=2):
    x=[1 if k in (i,j) else 0 for k in range(8)]
    for m in range(8):
        y=B(m)
        lhs=omul(f,omul(f,x,x),y); rhs=omul(f,x,omul(f,x,y))
        if lhs!=rhs: print('counterexample x = e%d+e%d, y=e%d:'%(i,j,m), lhs, rhs); break
    else: continue
    break
