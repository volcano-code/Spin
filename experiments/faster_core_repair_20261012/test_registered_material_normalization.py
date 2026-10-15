"""A nonzero exact face must not be deleted by rounded cross-product zero."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import unittest
import numpy as np
from make_it_stand import spin

HERE=Path(__file__).resolve().parent

def exact_chain(snapshot):
    points=snapshot['exact_material_vertices']
    chain=Counter()
    for face,birth in zip(snapshot['material_faces'],snapshot['material_birth_faces']):
        vertices=tuple(tuple(points[i]) for i in face)
        oriented=min(vertices[i:]+vertices[:i] for i in range(3))
        chain[(oriented,int(birth))]+=1
    return chain

class RegisteredMaterialNormalization(unittest.TestCase):
    def test_production_normalization_preserves_nonzero_exact_oriented_faces(self):
        registration=json.loads((HERE/'noodles-first-core-failures-registration-01/registration.json').read_text())
        config_bytes=Path(registration['mapping_config']).read_bytes()
        self.assertEqual(hashlib.sha256(config_bytes).hexdigest(),registration['mapping_config_sha256'])
        config=json.loads(config_bytes)
        payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        native,request=spin._native_request(spin.SpinItFasterRequest(**payload))
        settings=spin._native_evaluator_settings(native,spin.SpinEvaluatorConfig(**config['evaluator']))
        q=np.asarray(registration['failures']['material_boundary']['record']['coefficients'])
        result=native.capture_spin_it_faster_fullrank_geometry(request,q,settings)
        snapshots={s['stage']:s for s in result['snapshots']}
        before=snapshots['boolean-before-normalization'];after=snapshots['boolean-after-normalization']
        self.assertEqual(len(before['material_faces']),44856)
        self.assertEqual(len(after['material_faces']),44856,
                         'normalization removed a certified nonzero original face')
        self.assertEqual(exact_chain(before),exact_chain(after),
                         'same-call exact oriented face/source chain was changed')

if __name__=='__main__':unittest.main()
