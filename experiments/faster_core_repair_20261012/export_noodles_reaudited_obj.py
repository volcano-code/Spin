"""Export saved terminal arrays; retain explicit derived-geometry admission status."""
import hashlib
import json
from pathlib import Path
import numpy as np
from fractions import Fraction

HERE = Path(__file__).resolve().parent

def main():
    source = HERE/'noodles-terminal-recovery-job-03/output/material'
    receipt = json.loads((HERE/'noodles-terminal-bundle-reaudit-job-01/output/report.json').read_text())
    manifest = json.loads((source/'manifest.json').read_text())
    problem = json.loads((source/'problem.json').read_text())
    assert receipt['receipt']['status'] == 'MATERIAL_REAUDIT_IDENTICAL_FACTS'
    assert receipt['scientific_digest'] == manifest['scientific_digest']
    assert hashlib.sha256((source/'manifest.json').read_bytes()).hexdigest() == receipt['receipt']['source_manifest_sha256']
    out = HERE/'noodles-terminal-obj-01'
    out.mkdir(exist_ok=False)
    def load(ref):
        path = source/ref['name']
        assert path.parent == source and hashlib.sha256(path.read_bytes()).hexdigest() == ref['sha256']
        return np.load(path, allow_pickle=False)
    rows = {}
    meshes = {'original':problem['original'], 'fixed':problem['fixed_regions'][0]['mesh'],
              'design':problem['design_regions'][0]['mesh'], 'complete':manifest['final']['material'],
              'retained':manifest['final']['design_results'][0]['retained'],
              'cavity':manifest['final']['design_results'][0]['carved']}
    for name, mesh in meshes.items():
        v,f=load(mesh['vertices']),load(mesh['faces'])
        path=out/(name+'.obj')
        with path.open('x') as stream:
            stream.write('# metres; saved optimized material, no display transform\n')
            for p in v:stream.write('v '+' '.join(format(float(x),'.17g') for x in p)+'\n')
            for face in f:stream.write('f '+' '.join(str(int(x)+1) for x in face)+'\n')
        rv=[];rf=[]
        for line in path.read_text().splitlines():
            if line.startswith('v '):rv.append([float(x) for x in line.split()[1:]])
            elif line.startswith('f '):rf.append([int(x)-1 for x in line.split()[1:]])
        assert np.array_equal(rv,v) and np.array_equal(rf,f)
        cross=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
        candidates=np.flatnonzero(np.all(cross==0,axis=1))
        zero=[]
        for index in candidates:
            a,b,c=[[Fraction(float(x)) for x in v[j]] for j in f[index]]
            u=[b[i]-a[i] for i in range(3)];w=[c[i]-a[i] for i in range(3)]
            if all(u[i]*w[(i+1)%3]-u[(i+1)%3]*w[i]==0 for i in range(3)):zero.append(int(index))
        edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
        _, counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
        rows[name]=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            vertices=len(v),faces=len(f),roundtrip_identical=True,exact_zero_area_faces=zero,
            non_two_incident_edges=int(np.count_nonzero(counts!=2)),
            qualified_obj=False,remaining_checks=['self_intersection','bounded_authority_deviation','derived_partition'])
    report=dict(scientific_digest=manifest['scientific_digest'],length_unit=problem['length_unit'],
                mass_unit=problem['mass_unit'],physical=manifest['facts']['physical'],
                objects=rows,qualified_delivery=False)
    (out/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:{x:v[x] for x in ('faces','exact_zero_area_faces','non_two_incident_edges')} for k,v in rows.items()}))

if __name__=='__main__':main()
