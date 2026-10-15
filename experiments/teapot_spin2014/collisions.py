"""One shared teapot collision proxy pair, independent of carving strategy."""
import argparse
import json
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from teapot import write_mesh, write, file_hashes, check_files
from make_it_stand.pybullet_runner import Case,SimulationConfig,ensure_vhacd,VHACD_PRESETS


def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    a.output.mkdir(exist_ok=False)
    meta=json.loads((a.input/'manifest.json').read_text());check_files(meta['geometry_hashes'])
    scale=meta['meters_per_model_unit']; outer=trimesh.load(a.input/'outer.obj',force='mesh',process=False)
    outer.vertices*=scale; write_mesh(a.output/'outer-si.obj',outer)
    rod=trimesh.load(a.input/'rod.obj',force='mesh',process=False).convex_hull
    rod.vertices*=scale;write_mesh(a.output/'shaft-si.obj',rod)
    profiles={}
    for name,settings in VHACD_PRESETS.items():
        folder=a.output/name;folder.mkdir()
        print('V-HACD',name,flush=True)
        file,info=ensure_vhacd(Case('teapot',a.input/'teapot.obj',a.input/'teapot.obj'),folder,
            SimulationConfig(vhacd_resolution=settings['resolution'],vhacd_concavity=settings['concavity']))
        pieces=trimesh.load(file,force='mesh',process=False).split(only_watertight=False)
        paths=[a.output/'shaft-si.obj'];vertices=[rod.vertices]
        for i,part in enumerate(pieces):
            hull=part.convex_hull;hull.vertices*=scale
            path=folder/f'part-{i}.obj';write_mesh(path,hull);paths.append(path);vertices.append(hull.vertices)
        proxy=np.vstack(vertices);errors=[]
        for tilt in [0,2,5,15,45,90]:
            for direction in range(8):
                R=Rotation.from_rotvec(np.deg2rad(tilt)*np.array([np.cos(direction*np.pi/4),np.sin(direction*np.pi/4),0])).as_matrix()
                error=float((proxy@R.T)[:,2].min()-(outer.vertices@R.T)[:,2].min())
                errors.append(dict(tilt=tilt,direction=direction,error_m=error))
        profiles[name]=dict(files=[str(p.resolve()) for p in paths],hashes=file_hashes(paths),vhacd=info,
            support_errors=errors,near_upright_passed=max(abs(e['error_m']) for e in errors if e['tilt']<=5)<=.0001)
    write(a.output/'manifest.json',dict(outer=str((a.output/'outer-si.obj').resolve()),
        outer_hash=file_hashes([a.output/'outer-si.obj']),profiles=profiles,input_digest=meta['problem_digest'],
        source_hashes=file_hashes([Path(__file__),a.input/'manifest.json']),
        scope='SHAFT_CONVEX_BASE_PLUS_PADDED_VHACD_TEAPOT_LINKS_NOT_EXACT_CONTACT'))


if __name__=='__main__':main()
