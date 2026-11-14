"""Independent exact checker for C2--C3 using exterior minors.

Unlike verify_collision_lattices.py, this implementation constructs the
formal cubic resultant by a 3-by-3 permanent, tests equivalence on a
hyperplane by 2-by-2 minors, and tests equivalence on a line by
3-by-3 minors. Its integer-point test solves a 2-by-2 system over
all residue classes of a nonzero minor, not via a cross-product Bezout
base-point formula. Lines are identified by projective Pluecker
coordinates. This is independent of the primary checker routines.
"""
from collections import Counter
from itertools import combinations, permutations, product
from functools import reduce
from math import gcd

ZERO=(0,0,0,0)


def add4(a,b):return tuple(x+y for x,y in zip(a,b))
def sub4(a,b):return tuple(x-y for x,y in zip(a,b))
def shift4(a,k):return a[:3]+(a[3]+k,)

def primitive(v):
    g=reduce(gcd,v,0)
    assert g
    w=tuple(x//g for x in v)
    if next(x for x in w if x)<0:w=tuple(-x for x in w)
    return w


def polyadd_toggle(S,e):
    if e in S:S.remove(e)
    else:S.add(e)


def formal_cubic(residue):
    r,i,j=residue
    g=[set() for _ in range(3)]
    for ell, term in [(0,ZERO),(1,ZERO),(r,(1,0,0,0)),
                      (i,(0,1,0,0)),(j,(0,0,1,0))]:
        polyadd_toggle(g[ell],term)
    # Multiplication matrix from F_2[Y][X]/(X^3-Y)
    # [[g0,Yg2,Yg1],[g1,g0,Yg2],[g2,g1,g0]]
    mat=[[g[0],{shift4(e,1) for e in g[2]}, {shift4(e,1) for e in g[1]}],
         [g[1],g[0],{shift4(e,1) for e in g[2]}],
         [g[2],g[1],g[0]]]
    result=set()
    for rows in permutations(range(3)):
        entries=[mat[rows[col]][col] for col in range(3)]
        if any(not x for x in entries):continue
        for parts in product(*entries):
            e=add4(add4(parts[0],parts[1]),parts[2])
            polyadd_toggle(result,e)
    return sorted(result)


def plane_forms(forms):
    return sorted({primitive(sub4(a,b)) for a,b in combinations(forms,2)
                   if a[:3]!=b[:3]})


def minors2(a,b):
    return tuple(a[i]*b[j]-a[j]*b[i] for i,j in combinations(range(4),2))


def gw_hyperplane(forms,h):
    c=Counter(minors2(h,e) for e in forms)
    return sum(n%2 for n in c.values())


def det3(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def minors3(a,b,c):
    return tuple(det3([a[i] for i in cols], [b[i] for i in cols],
                      [c[i] for i in cols]) for cols in combinations(range(4),3))


def gw_line(forms,h1,h2):
    classes=Counter(minors3(h1,h2,e) for e in forms)
    return sum(n%2 for n in classes.values())


def cross3(a,b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def integer_line(h1,h2):
    a,b=h1[:3],h2[:3]
    nonzero=None
    for i,j in combinations(range(3),2):
        delta=a[i]*b[j]-a[j]*b[i]
        if delta:
            nonzero=(i,j,delta)
            break
    if nonzero is None:return None
    i,j,delta=nonzero
    free=3-i-j
    x0=None
    # Any integer solution has free-coordinate value in some class
    # modulo abs(delta); test every such class.
    for z in range(abs(delta)):
        rhs1=-h1[3]-a[free]*z
        rhs2=-h2[3]-b[free]*z
        num_i=rhs1*b[j]-rhs2*a[j]
        num_j=a[i]*rhs2-b[i]*rhs1
        if num_i%delta==0 and num_j%delta==0:
            x=[0]*3
            x[i]=num_i//delta;x[j]=num_j//delta;x[free]=z
            x0=tuple(x)
            break
    if x0 is None:return None
    v=primitive(cross3(a,b))
    p=x0+(1,)
    q=v+(0,)
    # Projective Pluecker coordinates of the affine line.
    key=primitive(minors2(p,q))
    return key,x0,v


def floor_div(a,b):
    return a//b

def ceil_div(a,b):
    return -((-a)//b)


def feasible(x0,v,residue):
    rr,ii,jj=residue
    T0,A0,B0=x0
    dT,dA,dB=v
    t0,t1=3*T0+rr,3*dT
    a0,a1=3*A0+ii,3*dA
    b0,b1=3*B0+jj,3*dB
    requirements=[(t1,t0-327),(b1,b0-2),
                  (a1-b1,a0-b0-1),(t1-a1,t0-a0-1)]
    lo,hi=None,None
    for slope,constant in requirements:
        if slope>0:
            lower=ceil_div(-constant,slope)
            lo=lower if lo is None else max(lo,lower)
        elif slope<0:
            upper=floor_div(-constant,slope)
            hi=upper if hi is None else min(hi,upper)
        elif constant<0:
            return None
    if lo is not None and hi is not None and lo>hi:return None
    parity=None
    if t1%2:
        parity=(1-t0)%2
    elif t0%2==0:
        return None
    if lo is not None:z=lo
    elif hi is not None:z=hi-2
    else:z=0
    if parity is not None and z%2!=parity:z+=1
    if hi is not None and z>hi:return None
    actual=(t0+t1*z,a0+a1*z,b0+b1*z)
    assert actual[0]>=327 and actual[0]%2 and 1<actual[2]<actual[1]<actual[0]
    return actual


def verify():
    formal_hist=Counter(); plane_hist=Counter()
    count_h=count_lines=count_low=0
    finalists=[]
    for residue in product(range(3),repeat=3):
        forms=formal_cubic(residue)
        formal_hist[len(forms)]+=1
        planes=plane_forms(forms)
        for h in planes:
            if gw_hyperplane(forms,h)<=5:
                count_h+=1
                plane_hist[h]+=1
        lines={}
        for h1,h2 in combinations(planes,2):
            got=integer_line(h1,h2)
            if got is None:continue
            key,x0,v=got
            # Generic weight uses minors, independent of base point.
            if key in lines:
                continue
            lines[key]=(x0,v,h1,h2)
        count_lines+=len(lines)
        for key,(x0,v,h1,h2) in lines.items():
            W=gw_line(forms,h1,h2)
            if W<=5:
                count_low+=1
                actual=feasible(x0,v,residue)
                if actual is not None:
                    finalists.append((residue,key,W,actual))
    assert formal_hist==Counter({13:18,14:7,17:2}),formal_hist
    assert count_h==75 and plane_hist==Counter({(0,0,1,0):16,(0,1,0,0):16,(1,0,0,0):16,(0,1,-1,0):9,(1,-1,0,0):9,(1,0,-1,0):9}),plane_hist
    assert (count_lines,count_low,len(finalists))==(11298,1026,7),(count_lines,count_low,len(finalists))
    assert all(x[2]==5 for x in finalists)
    expected={((0,0,0),(333,222,111)),((0,2,1),(327,326,325)),
              ((0,2,1),(327,164,163)),((0,2,2),(327,326,2)),
              ((1,1,1),(343,229,115)),((2,1,0),(329,4,3)),
              ((2,1,2),(329,247,83))}
    assert {(x[0],x[3]) for x in finalists}==expected,finalists
    print('Independent exterior-minor verifier C2/C3: PASS')
    print('formal support histogram',formal_hist)
    print('generic-weight<=5 hyperplanes',count_h,'by plane',plane_hist)
    print('integral affine lines',count_lines,'low-generic-weight',count_low,'admissible',len(finalists))
    for x in sorted(finalists):print('final',x[0],x[2],x[3])


if __name__=='__main__':
    verify()
