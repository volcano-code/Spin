import hashlib
import json
from pathlib import Path
import tempfile
import unittest


class SourceBindingTests(unittest.TestCase):
    def test_candidate_fixed_files_are_checked_before_construction(self):
        from clearance_source_binding import check_candidate_authority_inputs
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('heart-body.obj','pin.obj'):(root/name).write_text(name)
            frozen={name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                    for name in ('heart-body.obj','pin.obj')}
            frozen['spec.json']='spec-identity'
            identity={name:{'sha256':value} for name,value in frozen.items()}
            identity['spec']=identity.pop('spec.json')
            spec={'density':1200}
            construction={'frozen_outer_source':frozen,'inherited_physics_and_solver':spec}
            check_candidate_authority_inputs(root,construction,identity,spec)
            (root/'pin.obj').write_text('different fixed pin')
            with self.assertRaisesRegex(ValueError,'FIXED_GEOMETRY_CHANGED'):
                check_candidate_authority_inputs(root,construction,identity,spec)
            with self.assertRaisesRegex(ValueError,'INHERITED_SPEC'):
                check_candidate_authority_inputs(root,construction,identity,{'density':999})

    def test_rejects_source_spec_from_another_revision(self):
        from clearance_source_binding import check_source_contract
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'spec.json').write_text(json.dumps({'density_kg_m3':999.}))
            sha=hashlib.sha256((root/'spec.json').read_bytes()).hexdigest()
            (root/'preparation.json').write_text(json.dumps({'sources':{'spec.json':sha}}))
            with self.assertRaisesRegex(ValueError,'SOURCE_SPEC'):
                check_source_contract(root, {'density_kg_m3':1200.})

    def test_rejects_pin_bound_to_another_authority(self):
        from clearance_source_binding import check_fixed_source_identity
        source={'heart-body.obj':'body','pin.obj':'actual-pin','spec.json':'spec'}
        candidate=dict(source)
        recorded={name:{'sha256':sha} for name,sha in source.items()}
        recorded['pin.obj']['sha256']='different-pin'
        recorded['spec']=recorded.pop('spec.json')
        with self.assertRaisesRegex(ValueError,'AUTHORITY_SOURCE'):
            check_fixed_source_identity(source,candidate,recorded)


if __name__ == '__main__':unittest.main()
