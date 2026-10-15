"""Independent multiprecision geometry formulas and forward scalar dual algebra.

Development reference only. No production/native derivative or integrator import.
"""
import mpmath as mp

def classify(production,reference_low,reference_high,geometry_and_branch_valid):
    arrays=[list(mp.matrix(x)) for x in (production,reference_low,reference_high)]
    if not geometry_and_branch_valid or len({len(a) for a in arrays})!=1 or not arrays[0]:return 'INCONCLUSIVE'
    if not all(mp.isfinite(x) for a in arrays for x in a):return 'INCONCLUSIVE'
    if max(abs(a-b) for a,b in zip(arrays[1],arrays[2]))>mp.mpf('1e-20'):return 'INCONCLUSIVE'
    return 'SUPPORTED' if max(abs(a-b) for a,b in zip(arrays[0],arrays[2]))<=mp.mpf('1e-6') else 'INCONSISTENT'

def dot(a,b):return sum(x*y for x,y in zip(a,b))
def sub(a,b):return [x-y for x,y in zip(a,b)]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]

class Dual:
    def __init__(self,value,derivative=0):self.v=mp.mpf(value);self.d=mp.mpf(derivative)
    @staticmethod
    def cast(x):return x if isinstance(x,Dual) else Dual(x)
    def __add__(self,other):
        b=self.cast(other);return Dual(self.v+b.v,self.d+b.d)
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,-self.d)
    def __sub__(self,other):return self+-self.cast(other)
    def __rsub__(self,other):return self.cast(other)+-self
    def __mul__(self,other):
        b=self.cast(other);return Dual(self.v*b.v,self.d*b.v+self.v*b.d)
    __rmul__=__mul__
    def __truediv__(self,other):
        b=self.cast(other)
        if not b.v:raise ValueError('SINGULAR_DUAL_DIVISION')
        return Dual(self.v/b.v,(self.d*b.v-self.v*b.d)/(b.v*b.v))
    def __rtruediv__(self,other):return self.cast(other)/self

def quadratic(q):return mp.matrix([[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]])
def finite_point(q,normal):
    A=quadratic(q);b=mp.matrix(q[1:4]);n=mp.matrix(normal)
    u=mp.lu_solve(A,n);center=-mp.lu_solve(A,b)/2
    level=(b.T*mp.lu_solve(A,b))[0]/4-q[0]
    denominator=(n.T*u)[0]
    if level/denominator<=0:raise ValueError('NO_FINITE_POSITIVE_NORMAL_BRANCH')
    return list(center+mp.sqrt(level/denominator)*u)

def root_point_and_jacobian(q,edge_normal,radius,guess):
    def constraints(x,coefficients):
        A=quadratic(coefficients);b=mp.matrix(coefficients[1:4]);v=mp.matrix(x)
        gradient=2*A*v+b
        return [(v.T*A*v)[0]+dot(b,v)+coefficients[0],
                dot(v,v)-radius*radius,dot(edge_normal,gradient)]
    point=mp.findroot(lambda x,y,z:constraints([x,y,z],q),tuple(guess),
                      tol=mp.power(10,-mp.mp.dps+10),maxsteps=30)
    system=mp.matrix([[mp.diff(lambda t:constraints([point[k]+(t if k==j else 0) for k in range(3)],q)[i],0)
                      for j in range(3)] for i in range(3)])
    inverse=system**-1
    condition=mp.norm(system)*mp.norm(inverse)
    rhs=mp.matrix([[mp.diff(lambda t:constraints(point,[q[k]+(t if k==j else 0) for k in range(10)])[i],0)
                   for j in range(10)] for i in range(3)])
    return list(point),-inverse*rhs,condition

def line_plane(points):
    a,b,c,d,e=points
    n=cross(sub(d,c),sub(e,c));direction=sub(b,a)
    divisor=dot(n,direction)
    if not divisor.v:raise ValueError('SINGULAR_LINE_PLANE')
    t=dot(n,sub(c,a))/divisor
    return [x+t*delta for x,delta in zip(a,direction)]

def tetra_moments(vertices,faces,density):
    result=[Dual(0) for _ in range(10)]
    for face in faces:
        a,b,c=[vertices[i] for i in face]
        volume=dot(a,cross(b,c))*density/6
        total=[a[i]+b[i]+c[i] for i in range(3)]
        result[0]+=volume
        for i in range(3):result[i+1]+=volume*total[i]/4
        for k,(i,j) in enumerate(((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))):
            result[k+4]+=volume*(total[i]*total[j]+a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/20
    return result
