"""Measure an exact positive affine image in physical and normalized frames."""
from fractions import Fraction
import json
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import exact_fraction
from make_it_stand.verification.rational_text import fraction_from_rational_text,fraction_to_rational_text

HERE=Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    run=HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
    config=json.loads((run/'configuration.json').read_text())
    norm=json.loads((run/'result.json').read_text())['identity']['normalization']
    terminal=json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
    axes=[[Fraction(float(x)) for x in row] for row in config['request']['target_axes']]
    origin=[Fraction(float(x)) for x in norm['origin']]
    target_origin=[Fraction(float(x)) for x in config['request']['target_origin']]
    length=Fraction(float(norm['L']))
    determinant=(axes[0][0]*(axes[1][1]*axes[2][2]-axes[1][2]*axes[2][1])-
                 axes[0][1]*(axes[1][0]*axes[2][2]-axes[1][2]*axes[2][0])+
                 axes[0][2]*(axes[1][0]*axes[2][1]-axes[1][1]*axes[2][0]))
    if length<=0 or determinant<=0:raise ValueError('NOT_A_POSITIVE_INVERTIBLE_AFFINE_FRAME')
    rows=[];native=spin._native_module()
    for name,filename in (('initial','evaluation-000001.json'),('terminal',terminal)):
        material=json.loads((run/filename).read_text())
        physical=[[exact_fraction(x) for x in p] for p in material['material_proof']['exact_vertices']]
        target=[[fraction_from_rational_text(x) for x in p] for p in material['variable_rational_vertices']]
        replay=[]
        start=time.monotonic()
        for point in physical:
            model=[(point[i]-origin[i])/length-target_origin[i] for i in range(3)]
            replay.append([sum(model[i]*axes[i][j] for i in range(3)) for j in range(3)])
        affine_seconds=time.monotonic()-start
        if replay!=target:raise ValueError('SAME_CALL_PHYSICAL_TO_TARGET_MAP_NOT_EXACT')
        physical_tokens=[[fraction_to_rational_text(x) for x in p] for p in physical]
        target_tokens=[[fraction_to_rational_text(x) for x in p] for p in target]
        faces=np.asarray(material['variable_faces'],dtype=np.int64)
        start=time.monotonic();a=native.audit_exact_rational_triangle_boundary(physical_tokens,faces,'AnyNonzero')
        physical_seconds=time.monotonic()-start
        start=time.monotonic();b=native.audit_exact_rational_triangle_boundary(target_tokens,faces,'AnyNonzero')
        target_seconds=time.monotonic()-start
        expected=[fraction_from_rational_text(x)*determinant/length**3 for x in a['component_signed_volumes_exact']]
        actual=[fraction_from_rational_text(x) for x in b['component_signed_volumes_exact']]
        row=dict(name=name,faces=len(faces),affine_exact=True,affine_replay_seconds=affine_seconds,
            physical_audit_seconds=physical_seconds,target_audit_seconds=target_seconds,
            physical_valid=a['boundary_valid'],target_valid=b['boundary_valid'],
            component_signs_equal=a['component_volume_signs']==b['component_volume_signs'],
            exact_volume_scaling_equal=expected==actual,
            physical_pair_tests=a['exact_pair_tests'],target_pair_tests=b['exact_pair_tests'],
            maximum_physical_denominator_bits=max(x.denominator.bit_length() for p in physical for x in p),
            maximum_target_denominator_bits=max(x.denominator.bit_length() for p in target for x in p))
        rows.append(row);(output/'rows.json').write_text(json.dumps(rows))
    (output/'report.json').write_text(json.dumps(dict(rows=rows,
        passed=all(r['physical_valid'] and r['target_valid'] and r['component_signs_equal']
                   and r['exact_volume_scaling_equal'] for r in rows),
        determinant_exact=fraction_to_rational_text(determinant),
        volume_scale_exact=fraction_to_rational_text(determinant/length**3),
        boolean_calls=0,core_admission=False,
        scope='Fixed saved same-call affine geometry and performance diagnostic only')))


if __name__=='__main__':main(Path(sys.argv[1]))
