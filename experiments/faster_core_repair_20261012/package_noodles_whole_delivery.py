"""Publish only the independently checked whole/cavity representation."""
import json,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;out=HERE/'noodles-terminal-interior-patch-01'
read=lambda name:json.loads((out/name).read_text())
check=read('independent-output-check.json');chain=read('exact-surface-chain-check.json');partition=read('whole-partition-check.json')
assert check['boundary_valid'] and check['intersecting_pair_count']==0 and chain['passed']
assert partition['whole_cavity_outer_partition_passed']
rho=1200.;m=np.array(check['derived_moments'])*rho;mass=float(m[0]);center=m[1:4]/mass
xx,yy,zz,yz,xz,xy=m[4:];second=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
central=second-mass*np.outer(center,center);inertia=np.trace(central)*np.eye(3)-central
report=dict(schema='bounded-derived-whole-delivery-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',units='m',density_kg_m3=rho,mass_kg=mass,center_of_mass_m=center.tolist(),inertia_about_center_kg_m2=inertia.tolist(),files={name:dict(sha256=hashlib.sha256((out/name).read_bytes()).hexdigest()) for name in ['complete.obj','cavity.obj']},whole_boundary_geometry_passed=True,whole_cavity_outer_partition_passed=True,exact_surface_chain_equal=True,split_fixed_retained_partition_passed=False,complete_optimization_run_and_repeat_passed=False,physics_completed=False,manufacturing_certified=False,source_material='../noodles-terminal-recovery-job-03/output/material',note='Whole OBJ and complementary cavity only. Retained component candidate remains invalid; not a complete split assembly or full case acceptance.')
(out/'delivery.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
