from collections import Counter, defaultdict
from itertools import combinations,product
from math import gcd
from functools import reduce

def vec_add(a,b):return tuple(x+y for x,y in zip(a,b))
def vec_scale(c,a):return tuple(c*x for x in a)
def vec_sub(a,b):return tuple(x-y for x,y in zip(a,b))
def primitive(a):
 v=reduce(gcd,a,0)
 if not v: return None
 a=tuple(x//v for x in a)
 if next(x for x in a if x)!=abs(next(x for x in a if x)): a=tuple(-x for x in a)
 return a

def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))

def formal(triple):
 r,i,j=triple
 zero=(0,0,0,0)
 g=[[],[],[]]
 g[0].append(zero);g[1].append(zero)
 g[r].append((1,0,0,0));g[i].append((0,1,0,0));g[j].append((0,0,1,0))
 supp=Counter()
 for rr in range(3):
  for a in g[rr]:
   for b in g[rr]:
    c=list(vec_add(a,vec_scale(2,b)))
    c[-1]+=rr
    supp[tuple(c)]+=1
 for a in g[0]:
  for b in g[1]:
   for c in g[2]:
    l=list(vec_add(vec_add(a,b),c));l[-1]+=1;supp[tuple(l)]+=1
 return sorted(k for k,v in supp.items() if v%2)

def collision_hyperplanes(s):
 return sorted(set(primitive(vec_sub(a,b)) for a,b in combinations(s,2) if (a[0:3]!=b[0:3])))

def plane_gw(s,h):
 j=next(i for i in range(3) if h[i])
 reps=Counter(tuple(h[j]*v[k]-v[j]*h[k] for k in range(4)) for v in s)
 return sum(n%2 for n in reps.values())

def bezout(a,b):
 if b==0: return (abs(a),1 if a>0 else (-1 if a<0 else 0),0)
 z,x,y=bezout(b,a%b)
 return (z,y,x-(a//b)*y)

def bezout3(v):
 g,x,y=bezout(v[0],v[1])
 G,u,w=bezout(g,v[2])
 assert G==1
 res=(u*x,u*y,w)
 assert dot(v,res)==1
 return res

def candidate_line(h1,h2):
 d=cross(h1[:3],h2[:3]);g=reduce(gcd,d,0)
 if not g:return None
 m=tuple(h1[-1]*h2[i]-h2[-1]*h1[i] for i in range(3))
 if any(x%g for x in m):return None
 v=tuple(x//g for x in d);mo=tuple(x//g for x in m)
 if next(x for x in v if x)<0:
  v=tuple(-x for x in v);mo=tuple(-x for x in mo)
 # line key based only on v, m = x cross v
 key=(v,mo)
 return key

def base_point(v,mo):
 w=bezout3(v)
 x0=cross(w,mo)
 assert cross(x0,v)==mo
 return x0

def line_gw(s,v,x0):
 reps=Counter((dot(e[:3],v),dot(e[:3],x0)+e[3]) for e in s)
 return sum(x%2 for x in reps.values())

def feasible(x0,v,r):
 rr,ii,jj=r
 tn=3*v[0]; t0=3*x0[0]+rr
 an=3*v[1]; a0=3*x0[1]+ii
 bn=3*v[2]; b0=3*x0[2]+jj
 conditions=[(tn,t0-327),(bn,b0-2),(an-bn,a0-b0-1),(tn-an,t0-a0-1)]
 lo=None;hi=None
 for slope,const in conditions:
  if slope>0:
   bound= -((const)//slope) # ceil(-const/slope)
   lo=bound if lo is None else max(lo,bound)
  elif slope<0:
   bound= const//(-slope)
   hi=bound if hi is None else min(hi,bound)
  elif const<0:return False,None
 if lo is not None and hi is not None and lo>hi:return False,None
 if tn%2==0:
  if t0%2==0:return False,None
  parity=None
 else:parity=(1-t0)%2
 z=lo if lo is not None else ((hi-2) if hi is not None else 0)
 if parity is not None and z%2!=parity:z+=1
 if hi is not None and z>hi:return False,None
 return True,z

if __name__=='__main__':
 countH=0;histH=Counter(); total_lines=0;small_lines=0;final=[]
 formal_hist=Counter()
 for residue in product(range(3),repeat=3):
  s=formal(residue);formal_hist[len(s)]+=1
  H=collision_hyperplanes(s)
  for h in H:
   if plane_gw(s,h)<=5:
    countH+=1; histH[h]+=1
  lines=set()
  for h1,h2 in combinations(H,2):
   line=candidate_line(h1,h2)
   if line is not None:lines.add(line)
  total_lines+=len(lines)
  for v,mo in lines:
   x0=base_point(v,mo)
   W=line_gw(s,v,x0)
   if W<=5:
    small_lines+=1
    ok,z=feasible(x0,v,residue)
    if ok:
     T,A,B=tuple(x0[i]+z*v[i] for i in range(3))
     t,a,b=3*T+residue[0],3*A+residue[1],3*B+residue[2]
     assert t>=327 and t%2==1 and 1<b<a<t
     final.append((residue, v,mo, W, x0,z, (t,a,b)))
  print('residue',residue,'formal',len(s),'H',len(H),'lines',len(lines),flush=True)
 assert formal_hist==Counter({13:18,14:7,17:2}),formal_hist
 assert countH==75 and histH==Counter({(0,0,1,0):16,(0,1,0,0):16,(1,0,0,0):16,(0,1,-1,0):9,(1,-1,0,0):9,(1,0,-1,0):9}),histH
 assert (total_lines,small_lines,len(final))==(11298,1026,7),(total_lines,small_lines,len(final))
 assert all(x[3]==5 for x in final)
 expected = {((0,0,0),(3,2,1),(0,0,0)),((0,2,1),(1,1,1),(0,-1,1)),((0,2,1),(2,1,1),(0,-1,1)),((0,2,2),(1,1,0),(0,0,1)),((1,1,1),(3,2,1),(0,0,0)),((2,1,0),(1,0,0),(0,1,-1)),((2,1,2),(4,3,1),(1,-1,-1))}
 assert {(x[0],x[1],x[2]) for x in final}==expected
 print('C2/C3 exact assertions: PASS')
 print('FORMAL',formal_hist)
 print('C2',countH,'by-hyperplane',histH)
 print('C3',total_lines,small_lines,len(final))
 for x in sorted(final):print('SEVEN',x)
