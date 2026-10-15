"""Reference-derived preview only; never reads an optimized region or solver.

Recovers the visible top-layer outline under an affine camera approximation.
The central bore is assumed circular. Raster silhouettes contain projection,
antialiasing and extrusion-edge errors, so this is not original CAD recovery.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import cv2

def signed_area(loop):
    q=np.asarray(loop);return float(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1))/2)

def trace_outline(image, bearing_hint):
    mask=np.uint8(cv2.cvtColor(image,cv2.COLOR_BGR2GRAY)<175)*255
    contours,hierarchy=cv2.findContours(mask,cv2.RETR_CCOMP,cv2.CHAIN_APPROX_SIMPLE)
    if hierarchy is None:raise ValueError('NO_CONTOURS')
    hierarchy=hierarchy[0]
    outside=max((i for i in range(len(contours)) if hierarchy[i][3]<0),key=lambda i:cv2.contourArea(contours[i]))
    holes=[i for i in range(len(contours)) if hierarchy[i][3]==outside]
    if not holes:raise ValueError('NO_BEARING_OR_HOLES')
    bearing=min(holes,key=lambda i:np.linalg.norm(np.mean(contours[i][:,0],axis=0)-bearing_hint))
    center,diameters,angle=cv2.fitEllipse(contours[bearing])
    theta=np.deg2rad(angle);R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    # Symmetric inverse ellipse stretch preserves screen-right as model-right.
    affine=R@np.diag(2/np.array(diameters))@R.T
    ids=[outside,*holes];loops=[]
    for i in ids:
        points=cv2.approxPolyDP(contours[i],.5,True)[:,0].astype(float)
        q=(points-center)@affine.T;q[:,1]*=-1
        loops.append(q)
    scale=float(np.linalg.norm(loops[0],axis=1).max())
    loops=[q/scale for q in loops]
    for i,q in enumerate(loops):
        if (signed_area(q)>0)!=(i==0):loops[i]=q[::-1]
    return [q.tolist() for q in loops],dict(
        bearing_loop_index=ids.index(bearing),ellipse_center_px=list(center),
        ellipse_diameters_px=list(diameters),ellipse_angle_degrees=angle,
        image_to_plane_linear=(np.diag([1,-1])@affine/scale).tolist(),
        threshold=175,contour_simplification_px=.5,hole_count=len(holes),
        discarded_holes=0,postsolve_rotation=False,
        assumptions=['circular bearing bore','affine view approximation',
                     'visible silhouette approximates top face; thickness edges unresolved'])

def main(image_path,out):
    from solid_view import extrude
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MP
    from matplotlib.patches import PathPatch
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    image=cv2.imread(str(image_path))
    if image is None or image.shape[:2]!=(798,1418):raise ValueError('WRONG_REFERENCE_IMAGE')
    # Im9 is the UNOPTIMIZED Fig.9 rendering. Crop excludes the middle layer.
    source=image[:275];loops,record=trace_outline(source,(720,150))
    out.mkdir(parents=True,exist_ok=False)
    (out/'reference-unoptimized.png').write_bytes(image_path.read_bytes())
    record.update(classification='IMAGE_DERIVED_RECONSTRUCTION_PREVIEW',
        source_image_sha256=hashlib.sha256(image_path.read_bytes()).hexdigest(),
        source_pdf='reports/spin_faster_progress_20261011/build/sources/sif-final.pdf',
        page=10,image_name='Im9.jp2',crop_xywh=[0,0,1418,275],
        optimization_runs=0,optimization_domain_frozen=False,
        target_cut_position=None,meters_per_model_unit=.020,
        thickness='ASSUMED: bottom/top 0.02 units each; middle 0.12 units',
        loops=loops)
    meshes=[]
    for name,low,high in [('bottom',-.08,-.06),('middle',-.06,.06),('top',.06,.08)]:
        mesh=extrude(loops,low,high);mesh.export(out/(name+'-preview.obj'),digits=17)
        meshes.append(mesh)
    (out/'preview-input.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    from matplotlib import font_manager
    font_path='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
    font_manager.fontManager.addfont(font_path)
    plt.rcParams['font.family']=font_manager.FontProperties(fname=font_path).get_name()
    fig=plt.figure(figsize=(17,6),facecolor='#f3f5f7')
    ax=fig.add_subplot(131);ax.imshow(cv2.cvtColor(source,cv2.COLOR_BGR2RGB));ax.set_xlim(380,1060);ax.axis('off')
    ax.set_title('论文未优化上层：轮廓提取来源',fontsize=15)
    ax=fig.add_subplot(132);vertices=[];codes=[]
    for loop in loops:
        vertices.extend([*loop,loop[0]]);codes.extend([MP.MOVETO]+[MP.LINETO]*(len(loop)-1)+[MP.CLOSEPOLY])
    ax.add_patch(PathPatch(MP(vertices,codes),facecolor='#36434c',edgecolor='none'))
    ax.set(xlim=(-1.08,1.08),ylim=(-1.08,1.08));ax.set_aspect('equal');ax.axis('off');ax.set_title('恢复平面轮廓：保留全部可见孔洞',fontsize=15)
    # Use a depth buffer: Matplotlib's painter sorting can hide holes between
    # overlapping extruded layers. This is a geometry preview, not simulation.
    import pybullet as bullet
    client=bullet.connect(bullet.DIRECT)
    try:
        for name,offset,color in [('bottom',-.8,[.38,.43,.47,1]),('middle',0,[.81,.29,.57,1]),('top',.8,[.38,.43,.47,1])]:
            visual=bullet.createVisualShape(bullet.GEOM_MESH,fileName=str(out/(name+'-preview.obj')),rgbaColor=color)
            bullet.createMultiBody(baseMass=0,baseVisualShapeIndex=visual,baseCollisionShapeIndex=-1,basePosition=[0,0,offset])
        view=bullet.computeViewMatrix([.15,-5,3.2],[0,0,0],[0,0,1])
        projection=bullet.computeProjectionMatrixFOV(43,1,.01,20)
        pixels=bullet.getCameraImage(800,800,view,projection,renderer=bullet.ER_TINY_RENDERER,lightDirection=[-2,-3,5])[2]
    finally:bullet.disconnect(client)
    ax=fig.add_subplot(133);ax.imshow(np.asarray(pixels,np.uint8).reshape(800,800,4));ax.axis('off')
    ax.set_title('重建三层展开：尚未优化',fontsize=15)
    fig.text(.03,.055,'只读取未优化参考；圆孔去斜视是假设，厚度是重建参数。未预置空腔、未指定挖空象限。',fontsize=13,color='#36434c')
    fig.subplots_adjust(left=.02,right=.98,top=.88,bottom=.15,wspace=.05)
    fig.savefig(out/'preview.png',dpi=160);plt.close(fig)
    print(json.dumps({k:v for k,v in record.items() if k!='loops'},ensure_ascii=False))

if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
