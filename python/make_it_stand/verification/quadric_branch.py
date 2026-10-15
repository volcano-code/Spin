"""Independent rational invariants of a represented full-rank polynomial.

beta=a+b.x+x^T A x; order a,bx,by,bz,Axx,Ayy,Azz,Ayz,Axz,Axy.
No native classifier/eigenvectors, tolerance, mesh operation or solver call.
Analytic subtype equality does NOT prove discrete construction equivalence.
"""
from ._spin2014_authority import _arithmetic


def exact_full_rank_invariants(coefficients, *, arithmetic_backend='FRACTION_V1'):
    rational,read,write=_arithmetic(arithmetic_backend)
    if len(coefficients)!=10: raise ValueError('ten canonical rational coefficients required')
    q=[read(x) for x in coefficients]
    a=[[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]]
    work=[row[:] for row in a]
    positive=negative=0;k=0
    def swap(i,j):
        work[i],work[j]=work[j],work[i]
        for row in work: row[i],row[j]=row[j],row[i]
    while k<3:
        pivot=next((i for i in range(k,3) if work[i][i]!=0),None)
        if pivot is not None:
            swap(k,pivot);d=work[k][k]
            positive+=int(d>0);negative+=int(d<0)
            for i in range(k+1,3):
                for j in range(k+1,3): work[i][j]-=work[i][k]*work[k][j]/d
            k+=1
        else:
            pair=next(((i,j) for i in range(k,3) for j in range(i+1,3) if work[i][j]!=0),None)
            if pair is None: break
            i,j=pair;swap(k,i);swap(k+1,j)
            d=work[k][k+1];positive+=1;negative+=1
            for i in range(k+2,3):
                for j in range(k+2,3):
                    work[i][j]-=(work[i][k]*work[k+1][j]+work[i][k+1]*work[k][j])/d
            k+=2
    aa,b,c=a[0];_,d,e=a[1];_,_,f=a[2]
    adj=[[d*f-e*e,c*e-b*f,b*e-c*d],
         [c*e-b*f,aa*f-c*c,b*c-aa*e],
         [b*e-c*d,b*c-aa*e,aa*d-b*b]]
    determinant=aa*adj[0][0]+b*adj[0][1]+c*adj[0][2]
    result=dict(schema='independent-rational-quadric-invariants-v1',
        scope='ANALYTIC_POLYNOMIAL_NOT_DISCRETE_BRANCH_OR_OPTIMIZATION',
        signature=dict(positive=positive,negative=negative,zero=3-positive-negative),
        determinant=write(determinant),optimization_or_material_qualified=False)
    if determinant==0:
        return dict(result,centered_level=None,subtype='RANK_DEFICIENT_OUTSIDE_THIS_CHECK')
    level=q[0]-sum((q[i+1]*adj[i][j]*q[j+1] for i in range(3) for j in range(3)),rational(0))/(4*determinant)
    if level==0:
        subtype='CONE' if positive and negative else 'SINGLE_ZERO_POINT'
    elif positive==3 or negative==3:
        subtype='ELLIPSOID' if (positive==3 and level<0) or (negative==3 and level>0) else 'NO_REAL_ZERO_SET'
    else:
        matching=positive if level<0 else negative
        subtype='ONE_SHEET_HYPERBOLOID' if matching==2 else 'TWO_SHEET_HYPERBOLOID'
    return dict(result,centered_level=write(level),subtype=subtype)
