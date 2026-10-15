"""A trial's proof must be durably saved and selected by final state, not tail."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
from make_it_stand import spin
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'faster_heart_contact_20261011'))
from rounded_solver_journal import DurableFullRankMap


class DurableProof(unittest.TestCase):
    def test_one_real_mapping_writes_proof_and_binds_terminal_state(self):
        frozen=json.loads((Path(__file__).parent/'same-input-job-03/output/frozen-inputs.json').read_text())
        payload=dict(frozen['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
        request=spin.SpinItFasterRequest(**payload)
        settings=spin.SpinEvaluatorConfig(**dict(frozen['settings'],exact_binary64_material_boundary=True))
        with tempfile.TemporaryDirectory() as tmp:
            try: mapping=DurableFullRankMap(request,settings,Path(tmp)/'run',capture_material_proof=True)
            except TypeError as exc:self.fail('durable proof option unavailable: '+str(exc))
            with mapping:
                q=np.asarray(frozen['states']['ellipsoid'])
                value=mapping(q)
                saved=json.loads((mapping.directory/'evaluation-000001.json').read_text())
                proof=json.loads(saved['evaluation']['material_proof_json'])
                self.assertEqual(proof['geometry_hash'],value['geometry_hash'])
                # Synthetic terminal selectors exercise binding only, not solver convergence.
                result=dict(status='CONVERGED',quadric_coefficients=q.tolist(),trace=[dict(material_geometry_hash=value['geometry_hash'])])
                terminal=mapping.terminal_mapping(result)
                self.assertEqual(terminal['evaluation_index'],1)
                self.assertEqual(terminal['material_acceptance'],'NOT_CHECKED')
                bad=dict(result,quadric_coefficients=(q+0.01).tolist())
                with self.assertRaisesRegex(ValueError,'TERMINAL_MAPPING'):
                    mapping.terminal_mapping(bad)
                saved['evaluation']['geometry_hash']='tampered'
                (mapping.directory/'evaluation-000001.json').write_text(json.dumps(saved))
                with self.assertRaisesRegex(ValueError,'TERMINAL_MAPPING'):
                    mapping.terminal_mapping(result)


if __name__=='__main__': unittest.main()
