"""Freeze a new reconstructed flower and present two archived method cases."""
from pathlib import Path
import json,sys,hashlib,importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'multiple-cases-v1'
W=ROOT/'.superpowers/sdd/2026-10-07-spin2014-objectives-quadric'
sys.path[:0]=[str(HERE),str(ROOT/'python')]
from demo import authenticated_bytes

def dump(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False))

def freeze():
    OUT.mkdir(exist_ok=True)
    if any(OUT.iterdir()):raise FileExistsError('freeze output must be empty')
    t=np.arange(120)*2*np.pi/120;r=1+.20*np.cos(5*t)
    outer=np.column_stack([.12+r*np.cos(t),r*np.sin(t)])
    inner=np.column_stack([.12+.88*r*np.cos(t),.88*r*np.sin(t)])
    a=np.arange(64)*2*np.pi/64
    bearing=np.column_stack([.22*np.cos(a),.22*np.sin(a)])
    bore=np.column_stack([.12*np.cos(a),.12*np.sin(a)])
    # Radial ordered loops are simple and .88 homothetic about their star
    # center. A conservative inscribed-radius bound encloses the bearing.
    # The native polygon/layer validator remains the actual solver gate.
    assert .88*.8*np.cos(np.pi/120)-.12>.22
    problem=dict(design_region_loops=[inner.tolist(),bearing[::-1].tolist()],design_areal_density=.15,
        extrusion_interval=[-.075,.075],extrusion_volume_density=1.,
        fixed_extrusion_layers=[dict(loops=[outer.tolist(),inner[::-1].tolist()],z_interval=[-.075,.075],volume_density=1.),
                                dict(loops=[bearing.tolist(),bore[::-1].tolist()],z_interval=[-.075,.075],volume_density=1.)],
        maximum_iterations=100,maximum_function_evaluations=60,residual_tolerance=1e-10,
        initial_trust_radius=.25,maximum_trust_radius=1.)
    dump(OUT/'flower-cases.json',dict(schema='faster-planar-method-fixture-v1',
        scope='RECONSTRUCTED_FIVE_PETAL_METHOD_NOT_AUTHOR_FLOWER',
        protocol='Frozen before first solve; seeds0,1 input-bounds initialization; 8 attempts, 60 evaluations; no supplied solution',
        cases=[dict(id='reconstructed-five-petal-spinner-v1',mode='FIDGET_SPINNER',asset_class='RECONSTRUCTED',planar_problem=problem)]))
    dump(OUT/'flower-geometry.json',dict(outer=outer.tolist(),inner=inner.tolist(),bearing=bearing.tolist(),bore=bore.tolist(),
        label='RECONSTRUCTED_NOT_AUTHOR_ORIGINAL',meters_per_model_unit=.020,density_kg_m3=1200,geometry_frozen=True))
    print(OUT/'flower-cases.json')

def saved_case(name,path):
    import io,trimesh
    from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
    folder=path/'material-geometry';m=json.loads((folder/'manifest.json').read_text());meshes={};moments={};hashes={}
    for region,desc in m['meshes'].items():
        data={}
        for field in ('vertices','faces','rational_vertices'):
            record=desc.get(field)
            if record:
                payload=authenticated_bytes(folder,record);hashes[str(folder/record['name'])]=hashlib.sha256(payload).hexdigest()
                data[field]=json.loads(payload) if field=='rational_vertices' else np.load(io.BytesIO(payload),allow_pickle=False)
        meshes[region]=(data['vertices'],data['faces'])
        if region in ('outer','material'):
            moments[region]=[float(x) for x in integrate_exact_rational_raw_moments(data.get('rational_vertices',data['vertices'].tolist()),
                data['faces'].tolist(),arithmetic_backend='GMPY2_2_2_1_EXACT_V1')]
    np.testing.assert_allclose(moments['material'],m['identity']['audit']['total_moments_model'],atol=1e-12,rtol=1e-12)
    for region,(v,f) in meshes.items():
        # Display-only NPY remains separately identified; only outer collision
        # geometry is later used in dynamics.
        trimesh.Trimesh(v,f,process=False).export(OUT/f'{name}-{region}.obj',digits=17)
    dump(OUT/f'{name}-material.json',dict(source=str(path),classification='RECONSTRUCTED_METHOD_VALIDATION',
        moments=moments,source_sha256=hashes,manifest_identity=m['scientific_digest'],
        manufacturing_approved=False,inner_float_meshes='DISPLAY_ONLY'))
    return meshes,moments

def plot_saved():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib import font_manager
    from demo import to_physics
    plt.rcParams['font.family']=font_manager.FontProperties(fname='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc').get_name()
    plt.rcParams['axes.unicode_minus']=False
    rows=[('cuboid',W/'cuboid-yoyo-seed1-normalized-replay-run-01/output','Cuboid 悠悠球：非主轴输入 → 二次曲面挖空'),
          ('static',W/'static-lowering-normalized-replay-current1/output','静态平衡：平面是退化二次曲面')]
    for name,path,title in rows:
        meshes,moments=saved_case(name,path)
        fig=plt.figure(figsize=(13,5),facecolor='#edf2f8');fig.suptitle(title,fontsize=20)
        for k in (0,1):
            ax=fig.add_subplot(1,3,k+1,projection='3d',computed_zorder=False)
            v,f=meshes['outer'];ax.add_collection3d(Poly3DCollection(v[f],facecolors='#ced6e1' if k else '#d98b50',alpha=.16 if k else .65,zorder=1))
            if k:
                v,f=meshes['retained'];ax.add_collection3d(Poly3DCollection(v[f],facecolors='#149d98',alpha=1,zorder=3))
            outer=meshes['outer'][0];lo=outer.min(0);hi=outer.max(0)
            raw=np.array(moments['material' if k else 'outer']);c=raw[1:4]/raw[0]
            ax.scatter(*c,c='#dd3158',s=42,zorder=5)
            ax.plot([0,0],[0,0],[lo[2]-.15,hi[2]+.15],'--',color='#183e68',zorder=4)
            ax.set(xlim=(lo[0]-.1,hi[0]+.1),ylim=(lo[1]-.1,hi[1]+.1),zlim=(lo[2]-.1,hi[2]+.1),title='原始实体' if not k else '优化保留材料（外壳透明）')
            ax.set_box_aspect(hi-lo);ax.view_init(elev=23,azim=-62);ax.set_axis_off()
        ax=fig.add_subplot(1,3,3);ax.axis('off')
        before,after=[np.array(moments[x]) for x in ('outer','material')]
        c0=before[1:4]/before[0];c1=after[1:4]/after[0]
        if name=='cuboid':
            def misalignment(raw):
                I=to_physics(raw)[2];return np.linalg.norm(I[:2,2])/np.linalg.norm(I)
            text=f'目标：悠悠球惯量目标\n求解：Faster 二次曲面 / Newton\n\n质量：{before[0]:.4f} → {after[0]:.4f}\n非主轴残差：\n{misalignment(before):.3e} → {misalignment(after):.3e}\n\n实际两块设计域均保留\n未使用体素优化结果\n非作者原始 Cuboid 尺寸'
        else:
            h0=c0[2]+1.2;h1=c1[2]+1.2
            text=f'目标：2014 平方质心高度\n求解：Faster 平面专用路径\n\n支撑面 z = -1.2\n质心高度：{h0:.6f} → {h1:.6f}\n目标 x = 0.05\n实际 x：{c0[0]:.6f} → {c1[0]:.6f}\n\n重建多面体方法验证\n不是论文原始 tilted vase'
        ax.text(0,.9,text,va='top',fontsize=12,color='#142944')
        fig.text(.04,.035,'原始尺寸与参数未认证；重建案例用于展示论文功能。内部浮点网格仅显示。',fontsize=11)
        fig.subplots_adjust(left=.03,right=.98,top=.86,bottom=.09);fig.savefig(OUT/f'{name}.png',dpi=145);plt.close(fig)

if __name__=='__main__':
    if sys.argv[1]=='freeze':freeze()
    elif sys.argv[1]=='saved':plot_saved()
