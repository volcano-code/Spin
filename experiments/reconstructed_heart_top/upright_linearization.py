"""Local upright oracle for an axis-aligned (not necessarily symmetric) top."""
import numpy as np


def upright_matrix(mass,height,inertia,omega,gravity=9.81):
    I=np.asarray(inertia,float)
    if (I.shape!=(3,3) or not np.isfinite(I).all() or
        not np.isfinite([mass,height,omega,gravity]).all() or mass<=0 or height<=0 or gravity<0 or
        np.diag(I).min()<=0 or np.max(np.abs(I-np.diag(np.diag(I))))>1e-12*np.max(np.abs(I))):
        raise ValueError('requires positive diagonal pivot inertia and axial positive CoM height')
    i,j,k=np.diag(I); g=mass*gravity*height; w=omega
    # State: body components of world vertical nx,ny, transverse wx,wy.
    return np.array([[0,w,0,-1],[-w,0,1,0],[0,g/i,0,(j-k)*w/i],[-g/j,0,(k-i)*w/j,0]])
