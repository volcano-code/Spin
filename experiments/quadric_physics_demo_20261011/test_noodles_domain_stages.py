import unittest


class StageIdentityTests(unittest.TestCase):
    def test_stages_preserve_exact_frozen_exclusion_order(self):
        from noodles_domain_stages import next_batch
        names=['hub','pin',*['tube-'+str(i) for i in range(79)]]
        done=[]
        while len(done)<len(names):
            batch=next_batch(names,done,8)
            self.assertGreater(len(batch),0)
            self.assertLessEqual(len(batch),8)
            done+=batch
        self.assertEqual(done,names)

    def test_missing_or_reordered_predecessor_is_rejected(self):
        from noodles_domain_stages import next_batch
        for done in (['b'],['a','a'],['a','c']):
            with self.assertRaisesRegex(ValueError,'EXCLUSION_PREFIX_MISMATCH'):
                next_batch(['a','b','c'],done,8)

    def test_rational_checkpoint_preserves_actual_boolean_entity(self):
        import noodles_domain_stages
        import numpy as np
        import trimesh
        import json
        from pathlib import Path
        from tempfile import TemporaryDirectory
        import _make_it_stand as native
        from admit_rounded_input import exact_rows,construct,save_region
        def box(size,shift):
            mesh=trimesh.creation.box(extents=size);mesh.apply_translation(shift)
            return exact_rows(mesh.vertices),mesh.faces.tolist()
        domain=box([2,2,2],[0,0,0])
        a=box([1.2,1.1,3],[.2,.1,0]);b=box([.7,3,.8],[-.6,0,.1])
        first=construct(native,domain,a,'difference')
        expected=construct(native,first,b,'difference')
        with TemporaryDirectory() as directory:
            path=Path(directory)/'domain.obj';save_region(path,first)
            packet=json.loads(path.with_suffix('.rational.json').read_text())
            actual=construct(native,(packet['vertices'],packet['faces']),b,'difference')
        regions=[native.make_exact_material_region_from_rationals(v,np.asarray(f,dtype=np.int32)) for v,f in (expected,actual)]
        self.assertTrue(native.relate_exact_material_regions(*regions)['equal'])
