"""Synthetic regular tetrahedron: exercises the offline pipeline, not Heart."""
import hashlib
import json
from fractions import Fraction
import math
import pytest
import verify_heart_projective_reference as reference


def test_independent_pipeline_matches_analytic_tetrahedron_and_rejects_tampering(tmp_path,monkeypatch):
    monkeypatch.setattr(reference,'HERE',tmp_path)
    def save(name,data):
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(data));return path
    normals=[[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]]
    vertices=[[x/math.sqrt(3) for x in n] for n in normals]
    faces=[[0,2,1],[0,1,3],[0,3,2],[1,2,3]]
    volume=8/(9*math.sqrt(3))
    jacobian=[[0.]*10 for _ in range(10)]
    jacobian[0][0]=-1.5*volume
    for i in (4,5,6):jacobian[i][0]=-volume/6
    material=dict(quadric_coefficients=[-1,0,0,0,1,1,1,0,0,0],
        moment_jacobian=jacobian,variable_faces=faces,
        variable_rational_vertices=[[str(Fraction(x)) for x in row] for row in vertices],
        material_proof=dict(resolved_construction_records=[dict(vertex_index=i,point_type='EXPLICIT3D',
            operands=[dict(status='UNIQUE',source_vertices=[dict(input_object=1,input_vertex=i)])]) for i in range(4)]))
    saved=save('heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output/point-000-material.json',material)
    config=save('heart-cherchi-coordinate-map-job-01/output/registration.json',dict(
        length=1.,origin=[0.,0.,0.],original=dict(design_vertices=[]),
        normalized=dict(target_axes=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],target_origin=[0.,0.,0.],material_density=1.)))
    save('core-registered-failure-capture-job-01/output/source_side.json',dict(snapshots=[dict(source_fixed_normals=normals)]))
    save('heart-projective-reference-source-job-01/output/source.json',dict(
        paper_projective_hyperboloid=True,effective_icosphere_subdivisions=2,
        source_normal_vertices=[0,1,2,3],source_projective_rays=[[-1]*4 for _ in range(4)],
        interface_vertex_count=4,faces=faces,truncation_radius=3.,vertices=vertices))
    save('heart-projective-reference-source-job-01/output/report.json',dict(passed=True,
        material_sha256=hashlib.sha256(saved.read_bytes()).hexdigest(),
        registration_sha256=hashlib.sha256(config.read_bytes()).hexdigest()))
    reference.main(tmp_path/'result')
    report=json.loads((tmp_path/'result/report.json').read_text())
    assert report['status']=='SUPPORTED'
    assert report['wrong_jacobian_status']=='INCONSISTENT'
    assert not report['case_admission'] and not report['old_five_point_gate_passed']
    saved.write_text(saved.read_text()+' ')
    with pytest.raises(ValueError,match='IDENTITY_MISMATCH'):
        reference.main(tmp_path/'tampered')
