"""Derived CAD in the original model frame; not a new solver or material seal.

The input CAD sidecar is a local identity binding, not authentication. Actual
source/reloaded solid moments are checked against independent mesh integrals.
This does not certify exact set equality or manufacturing wall thickness.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from .mass_moments import TriangleMesh, compute_mass_properties
from .spin_common import canonical, digest, require


def _hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def _integral(mesh,density,origin=None,length=1.):
    v,f=mesh.arrays()
    if origin is not None:
        v=(v-origin)/length
    m=compute_mass_properties(TriangleMesh(v.tolist(),f.tolist()),density)
    q=m.second_moment
    return np.array([m.mass,*m.first_moment,q[0][0],q[1][1],q[2][2],q[1][2],q[0][2],q[0][1]])


def export_original_frame_cad(artifact,normalized_report,*,normalized_brep,
                              occt_executable,output_directory):
    artifact.revalidate()
    problem=artifact.problem; norm=problem.normalization
    unit=problem.length_unit
    require(unit in ('m','cm','mm'),'unsupported original CAD length unit')
    require(norm['axes']==np.eye(3).tolist(),'unsupported normalization axes')
    source=Path(normalized_brep).resolve(strict=True)
    executable=Path(occt_executable).resolve(strict=True)
    require(executable.is_file() and os.access(executable,os.X_OK),'invalid OCCT executable')
    require(normalized_report.get('status')=='FASTER_3D_MATERIAL_STEP_BREP_VALIDATED' and
            normalized_report.get('artifact_scientific_digest')==artifact.scientific_digest and
            normalized_report.get('artifact_evidence_digest')==artifact.evidence_digest and
            normalized_report.get('brep_sha256')==_hash(source),
            'normalized CAD identity does not bind the supplied material/file')
    retained={i:m for i,m,_ in artifact.design_results}
    regions=[(r.region_id,retained[r.region_id],r.density,'design') for r in problem.design_regions
             if not retained[r.region_id].explicit_empty]
    regions += [(r.region_id,r.mesh,r.density,'fixed') for r in problem.fixed_regions if not r.mesh.explicit_empty]
    require(regions and normalized_report.get('material_region_ids')==[r[0] for r in regions] and
            normalized_report.get('entity_roles')==[r[3] for r in regions],
            'CAD material labels/order differ from the sealed material')
    require(normalized_report.get('material_entity_ids')==
            [digest({'role':r[3],'region_id':r[0]}) for r in regions],
            'CAD entity identities differ from the sealed material')
    counts=normalized_report.get('entity_solid_counts')
    require(isinstance(counts,list) and len(counts)==len(regions) and
            all(type(c) is int and c>0 for c in counts),'invalid CAD solid grouping')
    require(sum(counts)==normalized_report.get('material_solid_count'),'CAD solid count mismatch')
    densities=np.array([r[2] for r in regions]); normalized_density=densities/norm['R']
    require(normalized_report.get('entity_densities')==normalized_density.tolist(),
            'CAD density differs from the sealed material')
    origin=np.array(norm['origin']); length=norm['L']
    expected=np.array([_integral(r[1],r[2]) for r in regions])
    expected_normalized=np.array([_integral(r[1],r[2]/norm['R'],origin,length) for r in regions])
    # Fixed input scales; neither coordinate offset nor the final solution can
    # enlarge this tolerance. All ten terms, including cross moments, checked.
    scales=norm['R']*length**np.array([3,*([4]*3),*([5]*6)])
    output=Path(output_directory);output.mkdir(parents=True,exist_ok=False)
    binding={'schema':'original-frame-cad-binding-v1','artifact_scientific_digest':artifact.scientific_digest,
        'normalized_cad_digest':digest(normalized_report),'source_brep_sha256':_hash(source),
        'occt_executable_sha256':_hash(executable),'normalization':norm,
        'length_unit':unit,'mass_unit':problem.mass_unit,'cad_frame':'ORIGINAL_MODEL',
        'material_region_ids':[r[0] for r in regions],'material_densities':densities.tolist(),
        'entity_solid_counts':counts,'scope':'DERIVED_CAD_MOMENT_VALIDATION_NOT_EXACT_SET_OR_MANUFACTURING'}
    (output/'request.json').write_text(canonical(binding)+'\n')
    stage=output/'staged';stage.mkdir()
    try:
        done=subprocess.run([str(executable),'--restore-material-frame',str(source),
            format(length,'.17g'),*[format(x,'.17g') for x in origin],unit,
            str(stage/'material.brep'),str(stage/'material.step')],
            capture_output=True,text=True,timeout=60)
        (output/'stdout.log').write_text(done.stdout);(output/'stderr.log').write_text(done.stderr)
        require(done.returncode==0,'original-frame CAD failed: '+done.stderr.strip())
        native=json.loads(done.stdout)
        require(isinstance(native,dict) and native.get('schema')=='original-cad-frame-v1' and native.get('length_unit')==unit,
                'unexpected original-frame CAD response')
        mm={'m':1000.,'cm':10.,'mm':1.}[unit]
        require(native.get('millimetres_per_unit')==mm,'CAD unit conversion mismatch')
        groups={};errors=[]
        for name in ('source','transformed','brep','step','step_mm'):
            solids=np.asarray(native['solid_geometric_moments'][name],dtype=float)
            require(solids.shape==(sum(counts),10) and np.isfinite(solids).all() and
                    (solids[:,0]>0).all(),'invalid reloaded CAD solid moments')
            pieces=[];cursor=0
            for count in counts:
                pieces.append(solids[cursor:cursor+count].sum(axis=0));cursor+=count
            measured=np.asarray(pieces)
            if name=='step_mm':
                measured/=mm**np.array([3,*([4]*3),*([5]*6)])
            measured*= (normalized_density if name=='source' else densities)[:,None]
            expected_values=expected_normalized if name=='source' else expected
            error=float(np.max(np.abs(measured-expected_values)/(np.ones(10) if name=='source' else scales)))
            require(error<=1e-6,'original-frame CAD moments differ from independent material: '+name)
            errors.append(error);groups[name]=measured.tolist()
        report={**binding,'entity_mass_moments':groups,
            'raw_mass_moments':{k:np.asarray(v).sum(axis=0).tolist() for k,v in groups.items()},
            'maximum_normalized_moment_error':max(errors),'native':native,
            'brep_sha256':_hash(stage/'material.brep'),'step_sha256':_hash(stage/'material.step'),
            'diagnostic_files':{name:{'sha256':_hash(output/name),'bytes':(output/name).stat().st_size}
                                for name in ('stdout.log','stderr.log')}}
        # Outputs are new and private; publish only after all reloads/integrals
        # passed. Preserve staged evidence on failure, never repair it silently.
        for name in ('material.brep','material.step'):
            os.link(stage/name,output/name)
        (output/'material.step.identity.json').write_text(canonical(report)+'\n')
        return report
    except (ValueError,RuntimeError,OSError,subprocess.SubprocessError,KeyError,TypeError) as error:
        (output/'failure.json').write_text(canonical({'type':type(error).__name__,'message':str(error)})+'\n')
        raise
