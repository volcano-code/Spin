"""Saved companion verifier; native whole-request checks remain mandatory."""
import json
from pathlib import Path
import gmpy2
from admit_bounded_request import text
from bounded_pin_correspondence import boundary_distance_certificate,numerical_pin_gate
from pin_role_witness import load_pin_boundary
from rational_request_view import sha


def load_companion(terminal_path,*,expected_source_hashes):
    terminal_path=Path(terminal_path);folder=terminal_path.parent
    record=json.loads(terminal_path.read_text())
    if record.get('schema')!='bounded-pin-companion-admission-v1' or \
            record.get('status')!='PIN_CORRESPONDENCE_CERTIFIED' or record.get('pin_correspondence_certified')is not True:
        raise ValueError('PIN_COMPANION_NOT_CERTIFIED')
    inputs=record.get('inputs',{});files=record.get('files',{})
    required=('pin','authority-outer','candidate-fixed.rational.json')
    if any(k not in inputs or k not in expected_source_hashes for k in required) or \
            any(k not in files for k in ('identity.json','derived-pin.rational.json')):
        raise ValueError('PIN_COMPANION_SOURCE_OR_FILES_REQUIRED')
    if any(key not in inputs or inputs[key]['sha256']!=wanted for key,wanted in expected_source_hashes.items()) or \
            any(sha(info['path'])!=info['sha256'] for info in inputs.values()):
        raise ValueError('PIN_COMPANION_SOURCE_IDENTITY_MISMATCH')
    for name,wanted in files.items():
        path=folder/name
        if path.parent.resolve()!=folder.resolve() or path.is_symlink() or sha(path)!=wanted:
            raise ValueError('PIN_COMPANION_FILE_IDENTITY_MISMATCH')
    pin=load_pin_boundary(inputs['pin']['path'])
    derived=json.loads((folder/'derived-pin.rational.json').read_text())
    if derived.get('role')!='BOUNDED_DERIVED_PIN_NOT_ORIGINAL_PIN' or \
            derived.get('construction')!='ORIGINAL_PIN_INTERSECT_SAME_DERIVED_FIXED' or \
            derived.get('original_pin_sha256')!=inputs['pin']['sha256'] or \
            derived.get('derived_fixed_sha256')!=inputs['candidate-fixed.rational.json']['sha256']:
        raise ValueError('PIN_COMPANION_ENTITY_ROLE_MISMATCH')
    outer=json.loads(Path(inputs['authority-outer']['path']).read_text())
    vertices=[[gmpy2.mpq(x) for x in row] for row in outer['vertices']]
    length=max(max(p[k] for p in vertices)-min(p[k] for p in vertices) for k in range(3))
    identity=json.loads((folder/'identity.json').read_text())
    if gmpy2.mpq(identity['normalization_length_authority_exact'])!=length:
        raise ValueError('PIN_COMPANION_AUTHORITY_LENGTH_MISMATCH')
    from rounded_input import exact_moments
    am=exact_moments(pin['vertices'],pin['faces']);dm=exact_moments(derived['vertices'],derived['faces'])
    boundary=boundary_distance_certificate(pin,derived,text(length/10**6))
    numerical=numerical_pin_gate(boundary,[text(x) for x in am],[text(x) for x in dm],
                                record['actual_removed_volume_exact'],text(length))
    if numerical['status']!='CERTIFIED':raise ValueError('PIN_COMPANION_RECOMPUTED_BOUND_NOT_CERTIFIED')
    return dict(derived_packet=derived,original_pin_packet=pin,recomputed_numerical_gate=numerical,
                attachment_path=str(terminal_path.resolve()),attachment_sha256=sha(terminal_path),
                native_fresh_checks_required=True,geometry_admission_granted=False,
                required_native_checks=['ORIGINAL_PIN_IN_AUTHORITY_FIXED','DERIVED_PIN_EQUALS_PIN_INTERSECT_SAME_FIXED',
                                        'DERIVED_PIN_IN_FIXED','DERIVED_PIN_IN_PIN','DERIVED_PIN_NO_DESIGN_OVERLAP'])
