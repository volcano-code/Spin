"""True connected heart/shaft union; rod is excluded from the design domain."""
from pathlib import Path
import subprocess,sys,tempfile
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/reconstructed_heart_top'))
from heart_top import heart_mesh


def write_mesh(path,mesh):
    # Binary64 roundtrip, unlike default 8-digit OBJ export.
    path.write_text(''.join('v '+' '.join(format(float(x),'.17g') for x in v)+'\n' for v in mesh.vertices)+
                    ''.join('f '+' '.join(str(int(x)+1) for x in f)+'\n' for f in mesh.faces))


def boolean(a,b,operation,work):
    with tempfile.TemporaryDirectory(prefix='cgal-',dir=work) as folder:
        folder=Path(folder); write_mesh(folder/'a.obj',a); write_mesh(folder/'b.obj',b)
        subprocess.run([str(ROOT/'build/heart-stem-contact-tools/mesh-boolean'),operation,
            str(folder/'a.obj'),str(folder/'b.obj'),str(folder/'out.obj')],check=True,timeout=180)
        if not (folder/'out.obj').stat().st_size:
            return trimesh.Trimesh(vertices=np.empty((0,3)),faces=np.empty((0,3),int),process=False)
        return trimesh.load(folder/'out.obj',force='mesh',process=False)


def rod_mesh():
    # Radius .025 model units (.625mm) at the rounded bottom; shaft radius
    # .065 (1.625mm). A finite triangulation, not an ideal mathematical point.
    angles=np.linspace(-np.pi/2,0,7)
    profile=[(0,-.16)]+[(.025*np.cos(t),-.135+.025*np.sin(t)) for t in angles[1:]]
    profile += [(.065,-.03),(.065,2.25),(0,2.25)]
    mesh=trimesh.creation.revolve(np.array(profile),sections=32)
    if mesh.volume<0: mesh.invert()
    return mesh


def make_input(work):
    heart=trimesh.Trimesh(*heart_mesh(),process=False); rod=rod_mesh()
    outer=boolean(heart,rod,'union',work)
    design0=heart.copy(); design0.vertices=[0,0,.95]+.72*(design0.vertices-[0,0,.95])
    design=boolean(design0,rod,'difference',work)
    fixed=boolean(outer,design,'difference',work)
    return dict(heart=heart,rod=rod,outer=outer,design=design,fixed=fixed)
