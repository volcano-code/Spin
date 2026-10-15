"""Render exact sealed meshes and actual simulated quaternion trajectories.

No geometry modification feeds back into optimization or validation. Cutaways
only hide triangles, and red depicts removed volume, not retained material.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial.transform import Rotation,Slerp


def setup(ax,title):
    ax.set(xlim=(-1.35,1.35),ylim=(-1.35,1.35),zlim=(-.2,2.5),title=title)
    ax.set_box_aspect((1,1,1)); ax.view_init(elev=12,azim=-78); ax.set_axis_off()
    ax.plot([0,0],[0,0],[0,2.2],color='#7c8399',lw=.8,ls='--')
    ax.scatter([0],[0],[0],color='#20242c',s=26)


def mesh(ax,v,f,color,alpha=1):
    t=v[f]; n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    n/=np.maximum(np.linalg.norm(n,axis=1)[:,None],1e-30)
    intensity=.60+.40*np.abs(n@np.array([.2,-.75,.63]))
    rgb=np.array(matplotlib.colors.to_rgb(color))
    obj=Poly3DCollection(t,facecolors=rgb*intensity[:,None],edgecolors='none',alpha=alpha)
    ax.add_collection3d(obj); return obj


def main():
    p=argparse.ArgumentParser(); p.add_argument('run',type=Path); p.add_argument('--video',action='store_true')
    a=p.parse_args(); dest=a.run/'visuals'; dest.mkdir(exist_ok=True)
    arrays=np.load(a.run/'validation/actual_meshes.npz')
    v,f=arrays['outer_v'],arrays['outer_f']; mv,mf=arrays['material_v'],arrays['material_f']
    cv,cf=arrays['cavity_v'],arrays['cavity_f']
    plt.rcParams.update({'font.size':11,'figure.facecolor':'#f6f5f1','axes.facecolor':'#f6f5f1'})
    fig=plt.figure(figsize=(13.5,5.5))
    titles=['Original solid heart','Final material + removed volume','Actual removed volume']
    for i,title in enumerate(titles):
        ax=fig.add_subplot(1,3,i+1,projection='3d'); setup(ax,title)
        if i==0: mesh(ax,v,f,'#6497d1')
        elif i==1:
            mesh(ax,mv,mf[(mv[mf].mean(axis=1)[:,0]<0)],'#e3a827')
            mesh(ax,cv,cf,'#dd5142')
        else: mesh(ax,cv,cf,'#dd5142')
    fig.suptitle('Self-built heart top | Actual Spin-It 2014 carving result',fontsize=17)
    fig.text(.5,.07,'Gold: retained material (display cutaway). Red: removed volume. Outer surface unchanged.',ha='center')
    fig.text(.5,.025,'RECONSTRUCTED example, not an original paper model. Mesh display is not simulation.',ha='center',color='#65676d')
    fig.subplots_adjust(left=0,right=1,bottom=.08,top=.82,wspace=0)
    fig.savefig(dest/'heart-carving.png',dpi=160); plt.close(fig)
    source=a.run/'supplemental-rotation'
    records={label:json.loads((source/f'{label}-35hz-dir0.json').read_text()) for label in ['original','carved']}
    fig,ax=plt.subplots(figsize=(10,4))
    for label,color in [('original','#537eb5'),('carved','#c68119')]:
        r=records[label]; ax.plot(r['time'],r['tilt_degrees'],label=label,color=color,lw=1)
    ax.set(xlabel='Actual simulation time (s)',ylabel='Axis tilt from vertical (degrees)',
           title='35 revolutions/s initial spin | 2 degree perturbation, direction 0')
    ax.legend(); ax.grid(alpha=.2); fig.text(.5,.01,'Ideal fixed pivot with gravity; no table contact or friction. Post-hoc speed selection.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1)); fig.savefig(dest/'actual-tilt-35hz.png',dpi=150); plt.close(fig)
    if not a.video: return
    fig=plt.figure(figsize=(10,5),dpi=100); objects=[]; interpolators=[]; labels=[]
    for i,(name,color) in enumerate([('original','#6497d1'),('carved','#e3a827')]):
        ax=fig.add_subplot(1,2,i+1,projection='3d'); setup(ax,name.capitalize())
        objects.append(mesh(ax,v,f,color))
        # Both display the unchanged outer surface; inertia comes from their
        # different material meshes. SLERP only interpolates solved orientations.
        r=records[name]; interpolators.append(Slerp(r['time'],Rotation.from_quat(r['quaternions'])))
        labels.append(ax.text2D(.15,.12,'',transform=ax.transAxes))
    title=fig.suptitle('Actual fixed-pivot dynamics | 35 revolutions/s',fontsize=15)
    fig.text(.5,.06,'0.1x slow motion | First 2 s of the 8 s run | Initial tilt 2 degrees',ha='center')
    fig.text(.5,.015,'Reconstructed heart. Ideal joint + gravity, no contact, no physical-print claim.',ha='center',fontsize=9)
    fig.subplots_adjust(left=0,right=1,bottom=.12,top=.82,wspace=0)
    writer=FFMpegWriter(fps=30,codec='libx264',extra_args=['-pix_fmt','yuv420p','-crf','22'])
    with writer.saving(fig,str(dest/'actual-heart-dynamics-35hz.mp4'),100):
        for j,t in enumerate(np.linspace(0,2,601)):
            for obj,slerp,label in zip(objects,interpolators,labels):
                R=slerp(t).as_matrix(); obj.set_verts((v@R.T)[f])
                tilt=np.rad2deg(np.arccos(np.clip(R[2,2],-1,1)))
                label.set_text(f't = {t:.3f} s   tilt = {tilt:.2f} deg')
            writer.grab_frame()
            if j%150==0: print(f'Actual trajectory frames {j}/600',flush=True)
    plt.close(fig)


if __name__=='__main__': main()
