"""Explicit compound collider: heart V-HACD plus exact-vertex shaft convex hull.

For a flat plane the convex hull has the same support function as its source
shaft. This does not assert identical contact manifolds or arbitrary obstacles.
"""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
import trimesh
from stem_geometry import ROOT,write_mesh
from run import write
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.pybullet_runner import Case,SimulationConfig,ensure_vhacd,VHACD_PRESETS


def main():
    p=argparse.ArgumentParser(); p.add_argument('run',type=Path)
    p.add_argument('--output-name',default='collision-tip-base'); a=p.parse_args()
    dest=a.run/a.output_name; dest.mkdir(exist_ok=False)
    config=json.loads((a.run/'config.json').read_text()); scale=config['meters_per_model_unit']
    original=trimesh.load(a.run/'outer.obj',force='mesh',process=False)
    outer=original.copy(); outer.vertices=np.asarray(outer.vertices)*scale
    write_mesh(dest/'outer-si.obj',outer)
    rod=trimesh.load(a.run/'rod.obj',force='mesh',process=False).convex_hull
    rod.vertices=np.asarray(rod.vertices)*scale; write_mesh(dest/'shaft-convex-si.obj',rod)
    records={}
    from scipy.spatial.transform import Rotation
    for preset,parameters in VHACD_PRESETS.items():
        folder=dest/preset; folder.mkdir()
        case=Case('heart',a.run/'heart.obj',a.run/'heart.obj')
        cfg=SimulationConfig(vhacd_resolution=parameters['resolution'],vhacd_concavity=parameters['concavity'])
        print('Decomposing heart',preset,flush=True)
        file,meta=ensure_vhacd(case,folder,cfg)
        combined=trimesh.load(file,force='mesh',process=False)
        # OBJ convex groups can share duplicate coordinates but no topology.
        parts=combined.split(only_watertight=False)
        # Only Bullet's base margin is adjustable in the installed build.
        paths=[str((dest/'shaft-convex-si.obj').resolve())]; vertices=[rod.vertices]
        for i,part in enumerate(parts):
            hull=part.convex_hull  # Explicitly match Bullet's per-part convexification.
            hull.vertices=np.asarray(hull.vertices)*scale
            path=folder/f'part-{i}.obj'; write_mesh(path,hull); paths.append(str(path.resolve()))
            vertices.append(hull.vertices)
        proxy=np.vstack(vertices); errors=[]
        for tilt in [0,2,5,15,45,90]:
            for azimuth in range(8):
                R=Rotation.from_rotvec(np.deg2rad(tilt)*np.array([np.cos(azimuth*np.pi/4),np.sin(azimuth*np.pi/4),0])).as_matrix()
                exact=float((outer.vertices@R.T)[:,2].min()); approx=float((proxy@R.T)[:,2].min())
                errors.append(dict(tilt=tilt,azimuth=azimuth,signed_support_error_m=approx-exact))
        max_error=max(abs(e['signed_support_error_m']) for e in errors if e['tilt']<=5)
        records[preset]=dict(files=paths,mesh_sha256={f:hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in paths},
            heart_vhacd=meta,support_audit=errors,near_upright_max_support_error_m=max_error,
            near_upright_support_passed=max_error<=.0001)
    write(dest/'manifest.json',dict(scope='COMPOUND_HEART_VHACD_PLUS_SHAFT_VERTEX_CONVEX_HULL',
        base_part='shaft',heart_link_margin_policy='MEASURE_RUNTIME_PADDING_NOT_EXACT_CONTACT',
        source_outer_sha256=hashlib.sha256((a.run/'outer.obj').read_bytes()).hexdigest(),
        meters_per_model_unit=scale,
        outer_file=str((dest/'outer-si.obj').resolve()),outer_sha256=hashlib.sha256((dest/'outer-si.obj').read_bytes()).hexdigest(),
        profiles=records,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print('Collision preparation finished',flush=True)


if __name__=='__main__': main()
