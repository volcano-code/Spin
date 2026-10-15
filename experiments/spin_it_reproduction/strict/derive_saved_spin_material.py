#!/usr/bin/env python3
"""Reaudit a bound saved material and export a NEW checked float derivative.

No optimizer, CAD repair or manufacturing approval. Geometry is authenticated
and parsed by the existing bundle reader; historical solver freshness is never
upgraded. Use the job supervisor for resource limits on large saved objects.
"""
import argparse
import os
from pathlib import Path

from make_it_stand.verification import derive_float_material
from make_it_stand.verification.spin_bundle import reaudit_spin_artifact, save_spin_artifact
from make_it_stand.verification.spin_common import canonical, require
from make_it_stand.verification.run_evidence import _path, _record


def run(source,output,*,expected_problem_digest,expected_scientific_digest):
    root=Path(output).absolute()
    _path(root.parent,directory=True)
    if root.exists() or root.is_symlink():
        raise FileExistsError(str(root))
    # Producer observations are separate from geometry identity and are not
    # compilation or loaded-memory attestation. The material gate is unchanged.
    from make_it_stand.verification import derived_material, _coplanar_boundary
    producers=[Path(__file__),Path(derived_material.__file__),Path(_coplanar_boundary.__file__)]
    before={p.name:_record(p,label=p.name) for p in producers}
    root.mkdir(mode=0o700)
    authority,reaudit=reaudit_spin_artifact(source,
        expected_problem_digest=expected_problem_digest,
        expected_scientific_digest=expected_scientific_digest)
    candidate,relation=derive_float_material(authority)
    saved=save_spin_artifact(candidate,root/'material')
    require(saved.scientific_digest==relation['derived_scientific_digest'],
            'saved material differs from derived relation')
    after={p.name:_record(p,label=p.name) for p in producers}
    require(before==after,'output producer changed during derivation')
    record={'schema':'saved-derived-spin-material-v1',
        'scope':'SAVED_MATERIAL_TRANSFORMATION_NOT_SOLVER_OR_MANUFACTURING_APPROVAL',
        'source_reaudit':reaudit,'relation':relation,
        'producer_observation':{'scope':'OBSERVED_FILES_ON_DISK_ONLY','files':after},
        'material_manifest':_record(root/'material/manifest.json',label='material/manifest.json')}
    # The receipt is last and exclusive. Interrupted/failed exports retain all
    # partial payloads but cannot leave a completed derivation receipt.
    with (root/'derivation.json').open('x',encoding='ascii') as stream:
        stream.write(canonical(record)+'\n');stream.flush();os.fsync(stream.fileno())
    # Persist both the receipt/material entries and the newly created root
    # entry in its parent before acknowledging a complete output.
    for directory in (root,root.parent):
        fd=os.open(directory,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW)
        try:os.fsync(fd)
        finally:os.close(fd)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--expected-problem-digest',required=True)
    parser.add_argument('--expected-scientific-digest',required=True)
    args=parser.parse_args()
    record=run(args.source,args.output,expected_problem_digest=args.expected_problem_digest,
               expected_scientific_digest=args.expected_scientific_digest)
    print(canonical({'derived_scientific_digest':record['relation']['derived_scientific_digest'],
                     'manufacturing_admission':record['relation']['manufacturing_admission'],
                     'solver_freshness':record['source_reaudit']['solver_freshness']}))


if __name__=='__main__':main()
