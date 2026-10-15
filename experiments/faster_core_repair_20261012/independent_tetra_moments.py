"""Independent signed tetrahedra, using exact GMP rational arithmetic."""
import gmpy2
import numpy as np


def integrate(vertices,faces,density):
    points=[]
    for row in vertices:
        points.append([gmpy2.mpq(x.numerator,x.denominator) if hasattr(x,'numerator') else gmpy2.mpq(x) for x in row])
    result=[gmpy2.mpq(0) for _ in range(10)]
    rho=gmpy2.mpq(float(density))
    for triangle in faces:
        a,b,c=[points[int(index)] for index in triangle]
        volume=(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6
        total=[a[i]+b[i]+c[i] for i in range(3)]
        result[0]+=rho*volume
        for i in range(3):result[i+1]+=rho*volume*total[i]/4
        for index,(i,j) in enumerate(((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))):
            result[index+4]+=rho*volume*(total[i]*total[j]+a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/20
    return np.asarray([float(x) for x in result])
