"""Self-built heart and independent ideal-pivot dynamics. Not paper-original.

Geometry is in model units; runner declares the fixed SI conversion. No contact,
friction, aerodynamic drag, or physical-print claim is made by this simulator.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation


def heart_mesh(samples=48, layers=12):
    if samples < 16 or samples % 4 or layers < 4 or layers % 2:
        raise ValueError('need even layers and at least16 profile samples divisible by4')
    t = np.arange(samples)*2*np.pi/samples
    x = 16*np.sin(t)**3/17
    z = (13*np.cos(t)-5*np.cos(2*t)-2*np.cos(3*t)-np.cos(4*t))/17+1
    x[np.abs(x)<1e-14]=0; z[np.abs(z)<1e-14]=0
    v=[]; faces=[]
    for a in np.linspace(-np.pi/2,np.pi/2,layers+1)[1:-1]:
        v.extend(zip(np.cos(a)*x,np.full(samples,.38*np.sin(a)),.9+np.cos(a)*(z-.9)))
    for r in range(layers-2):
        for j in range(samples):
            a=r*samples+j; b=r*samples+(j+1)%samples; c=a+samples; d=b+samples
            faces.extend(((a,b,d),(a,d,c)))
    lower=len(v); v.append((0,-.38,.9)); upper=len(v); v.append((0,.38,.9))
    for j in range(samples):
        k=(j+1)%samples; last=(layers-2)*samples
        faces.extend(((lower,k,j),(upper,last+j,last+k)))
    v=np.asarray(v,dtype=float); f=np.asarray(faces,dtype=np.int32)
    if integrate(v,f)[0]<0: f=f[:,::-1].copy()
    return v,f


def integrate(v,f):
    """Signed tetrahedral volume, first and full second geometric moments."""
    t=np.asarray(v)[np.asarray(f)]; a,b,c=t[:,0],t[:,1],t[:,2]
    volume=np.einsum('ij,ij->i',a,np.cross(b,c))/6
    sums=t.sum(axis=1)
    first=np.einsum('i,ij->j',volume,sums)/4
    second=np.einsum('i,ijk->jk',volume,
        np.einsum('ni,nj->nij',sums,sums)+np.einsum('nki,nkj->nij',t,t))/20
    return np.array([volume.sum(),*first,*np.diag(second),second[1,2],second[0,2],second[0,1]])


def properties(s):
    s=np.asarray(s,float)
    if s.shape!=(10,) or not np.isfinite(s).all() or s[0]<=0:
        raise ValueError('finite positive material mass required')
    q=np.array([[s[4],s[9],s[8]],[s[9],s[5],s[7]],[s[8],s[7],s[6]]])
    return s[0],s[1:4]/s[0],np.trace(q)*np.eye(3)-q


def simulate(mass,center,inertia,omega,*,duration=8.,gravity=9.81,
             tilt_degrees=2.,azimuth=0.,max_step=1/240):
    """Euler equations about ideal fixed pivot + body-to-world quaternion.

    The pivot supplies force, no torque. Gravity acts at actual material CoM.
    Energy and vertical angular momentum should be conserved; full world angular
    momentum is an oracle only for gravity=0. All units are SI here.
    """
    c=np.asarray(center,float); I=np.asarray(inertia,float); w=np.asarray(omega,float)
    if (c.shape!=(3,) or w.shape!=(3,) or I.shape!=(3,3) or
        not np.isfinite(np.r_[mass,c,w,I.ravel(),duration,gravity,max_step,tilt_degrees,azimuth]).all()
        or mass<=0 or duration<=0 or max_step<=0 or gravity<0
        or not np.allclose(I,I.T,atol=1e-15) or np.linalg.eigvalsh(I).min()<=0):
        raise ValueError('invalid physical state/inertia')
    inv=np.linalg.inv(I)
    q=Rotation.from_rotvec(np.deg2rad(tilt_degrees)*np.array([np.cos(azimuth),np.sin(azimuth),0])).as_quat()
    def rhs(_,state):
        q,w=state[:4],state[4:]
        R=Rotation.from_quat(q).as_matrix()
        torque=np.cross(c,R.T@np.array([0.,0.,-mass*gravity]))
        dq=.5*np.r_[q[3]*w+np.cross(q[:3],w),-q[:3]@w]
        return np.r_[dq,inv@(torque-np.cross(w,I@w))]
    times=np.linspace(0,duration,int(np.ceil(duration*240))+1)
    sol=solve_ivp(rhs,(0,duration),np.r_[q,w],method='DOP853',rtol=1e-10,atol=1e-12,
                  max_step=max_step,t_eval=times)
    if not sol.success or not np.isfinite(sol.y).all(): raise RuntimeError(sol.message)
    q=sol.y[:4].T; w=sol.y[4:].T; R=Rotation.from_quat(q).as_matrix()
    momentum=np.einsum('nij,nj->ni',R,w@I.T)
    potential=mass*gravity*np.einsum('nij,j->ni',R,c)[:,2]
    kinetic=.5*np.einsum('ni,ni->n',w,w@I.T); energy=kinetic+potential
    tilt=np.rad2deg(np.arccos(np.clip(R[:,2,2],-1,1)))
    e_scale=max(abs(energy[0]),abs(kinetic[0]),1e-30)
    l_scale=max(np.linalg.norm(momentum[0]),1e-30)
    return dict(time=times.tolist(),quaternions=q.tolist(),angular_velocity_body=w.tolist(),
        tilt_degrees=tilt.tolist(),
        max_tilt_degrees=float(tilt.max()),final_tilt_degrees=float(tilt[-1]),
        relative_energy_drift=float(np.max(np.abs(energy-energy[0]))/e_scale),
        relative_momentum_drift=float(np.max(np.linalg.norm(momentum-momentum[0],axis=1))/l_scale),
        relative_vertical_momentum_drift=float(np.max(np.abs(momentum[:,2]-momentum[0,2]))/l_scale),
        function_evaluations=sol.nfev)
