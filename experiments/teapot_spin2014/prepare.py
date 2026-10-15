"""Lock a genuine Utah mesh, explicitly close its four parts, and form solids."""
import argparse
from pathlib import Path
import time
import traceback
import numpy as np
import trimesh
from teapot import (ROOT, SOURCE_SHA256, boolean, cap_planar_boundaries, canonical_indices, shaft_mesh,
                    write_mesh, write, file_hashes, check_files)


def main():
    p = argparse.ArgumentParser(); p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    out = a.output.resolve(); out.mkdir(exist_ok=False); start = time.monotonic()
    sources = file_hashes([a.source, Path(__file__), Path(__file__).with_name('teapot.py'),
        ROOT/'experiments/heart_stem_contact/stem_geometry.py', ROOT/'build/heart-stem-contact-tools/mesh-boolean'])
    write(out/'source.json', dict(files=sources))
    try:
        if sources[str(a.source.resolve())] != SOURCE_SHA256: raise ValueError('UNRECOGNIZED_SOURCE')
        raw = trimesh.load(a.source, force='mesh', process=False)
        # Exact coordinate deduplication only; never round, smooth or remesh.
        v, inverse = np.unique(raw.vertices, axis=0, return_inverse=True)
        merged = trimesh.Trimesh(v, inverse[raw.faces], process=False)
        parts = sorted(merged.split(only_watertight=False), key=lambda x: -len(x.faces))
        if len(parts) != 4: raise ValueError('SOURCE_COMPONENT_COUNT_CHANGED')
        closed = []; audits = []
        for i, part in enumerate(parts):
            mesh, audit = cap_planar_boundaries(part)
            # y-up source -> z-up, proper rigid rotation then uniform scale.
            xyz = mesh.vertices.copy()
            mesh.vertices = .4*xyz[:,[0,2,1]]*np.array([1,-1,1])+[0,0,.15]
            write_mesh(out/f'closed-part-{i}.obj', mesh)
            closed.append(mesh); audits.append(audit)
        pot = closed[0]
        for part in closed[1:]: pot = boolean(pot, part, 'union', out)
        rod = shaft_mesh(); outer = boolean(pot, rod, 'union', out)
        design = closed[0].copy()
        design.vertices = [0,0,.67]+.84*(design.vertices-[0,0,.67])
        design = boolean(design, rod, 'difference', out)
        design, index_audit = canonical_indices(design)
        fixed = boolean(outer, design, 'difference', out)
        for name, mesh in [('teapot',pot),('rod',rod),('outer',outer),('design',design),('fixed',fixed)]:
            write_mesh(out/(name+'.obj'), mesh)
        from make_it_stand.verification.spin_problem import SpinArtifactProblem
        arrays = lambda m: (np.asarray(m.vertices),np.asarray(m.faces))
        problem = SpinArtifactProblem(original=arrays(outer), fixed_regions=[('fixed-shell-and-shaft',arrays(fixed),1.)],
            design_regions=[('body-interior',arrays(design),1.)],length_unit='25mm',mass_unit='18.75g',
            target={'mode':'SPINNING_TOP','axis':[0,0,1],'origin':[0,0,0],'support_points':[]}, minimum_outer_wall=.015)
        write(out/'input.json',problem.identity)
        check_files(sources)
        write(out/'manifest.json',dict(status='LOCKED_RECONSTRUCTED_INPUT',original_paper_case=False,
            source_url='https://graphics.stanford.edu/courses/cs148-10-summer/as3/code/as3/teapot.obj',
            asset_evidence='SOURCE_MODEL_ONLY_THEN_EXPLICIT_RECONSTRUCTION',parameter_evidence='IMPLEMENTATION_DEFAULT',
            source_vertices=len(raw.vertices),source_triangles=len(raw.faces),exact_unique_vertices=len(v),caps=audits,
            design_index_canonicalization=index_audit,
            operations=['exact coordinate weld','planar caps','proper axis rotation','uniform scale .4',
                        'CGAL component union','rounded shaft union','body scale .84 design minus shaft'],
            meters_per_model_unit=.025,density_kg_m3=1200.,minimum_outer_wall=.015,
            geometry_hashes=file_hashes([out/(name+'.obj') for name in ['outer','design','fixed','rod','teapot']]),
            problem_digest=problem.input_digest,source_files=sources,elapsed_seconds=time.monotonic()-start))
        print('LOCKED',problem.input_digest,flush=True)
    except Exception as e:
        write(out/'failure.json',dict(error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        raise


if __name__ == '__main__': main()
