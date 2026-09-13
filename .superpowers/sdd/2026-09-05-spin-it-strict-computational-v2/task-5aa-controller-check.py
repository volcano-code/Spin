"""Independent all-component coverage and exact distance checks on raw outputs.
No native imports; fixed synthetic fixtures only, no seal or paper authority.
"""
from collections import defaultdict
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path
import sys

spec=importlib.util.spec_from_file_location('region_oracle',Path(__file__).with_name('task-5z-controller-check.py'))
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
binary=prior.binary;geo=prior.geo
IDS={'empty','one_cavity','two_cavities','two_cavities_rejected','threshold_rejected','resource_limited','not_admitted'}


def connected_faces(mesh):
    edge_faces=defaultdict(list)
    for i,face in enumerate(mesh['faces']):
        for a,b in zip(face,face[1:]+face[:1]):edge_faces[tuple(sorted((a,b)))].append(i)
    adjacency=defaultdict(set)
    for pair in edge_faces.values():
        assert len(pair)==2
        a,b=pair;adjacency[a].add(b);adjacency[b].add(a)
    unseen=set(range(len(mesh['faces'])));groups=[]
    while unseen:
        todo=[min(unseen)];visited=set()
        while todo:
            i=todo.pop()
            if i in visited:continue
            visited.add(i);todo.extend(adjacency[i]-visited)
        unseen-=visited;groups.append(sorted(visited))
    return groups


def verify(row):
    assert row['id'] in IDS and row['scope']=='KNOWN_SOLUTION_VALIDATION'
    assert row['sealed'] is False and row['paper_reproduction'] is False
    assert row['metric_scope']=='BOUNDARY_COMPONENT_SEPARATION' and row['full_thickness']=='NOT_ESTABLISHED'
    assert row['input_hash']==row['upstream_hash'] and len(row['scientific_hash'])==64
    assert row['minimum']==(.2 if row['id'] in {'threshold_rejected','two_cavities_rejected'} else .05)
    assert row['nef_limit']==(1 if row['id']=='resource_limited' else 1000000)
    if row['id']=='not_admitted':
        assert row['upstream_validated'] is False and row['admitted'] is False
        assert row['outcome']==1 and row['published'] is False and row['same_input'] is False
        return 0.
    assert row['upstream_validated'] is True and row['admitted'] is True
    if row['id']=='resource_limited':
        assert row['outcome']==3 and row['published'] is False and row['same_binary'] is False
        assert row['nef_exceeded'] is True
        return 0.
    assert row['identities'] is True and row['coverage'] is True
    material,s=row['material'],row['source']
    error=binary.source.check_source(s)
    assert s['handles']==26 and s['density']==1.17
    expected=binary.boundary_chain(s['points'],s['outer_faces'])
    expected.update(binary.oriented_triangle(t,(0,2,1)) for t in binary.selected_boundary(s).elements())
    assert binary.boundary_chain(material['vertices'],material['faces'])==expected
    world=binary.integrate(material,s['density'])
    assert max(abs(world[i]-s['moments'][i]) for i in range(10))/s['density'] <= 1e-10
    assert row['material_input']==prior.mesh_digest(material)
    groups=connected_faces(material)
    count=1 if row['id']=='empty' else 3 if row['id'].startswith('two_cavities') else 2
    assert row['component_count']==len(groups)==len(row['components'])==count
    labels=[None]*len(material['faces']);meshes=[]
    for i,(faces,c) in enumerate(zip(groups,row['components'])):
        assert c['id']==i and c['faces']==faces
        for j in faces:labels[j]=i
        vertices=sorted({j for f in faces for j in material['faces'][f]})
        assert c['vertices']==vertices
        remap={j:k for k,j in enumerate(vertices)}
        mesh={'vertices':[material['vertices'][j] for j in vertices],
              'faces':[[remap[j] for j in material['faces'][f]] for f in faces]}
        raw=geo.moments([geo.vec(v) for v in mesh['vertices']],mesh['faces'])
        sign=1 if raw[0]>0 else -1
        assert raw[0]!=0 and c['sign']==sign== (1 if i==0 else -1)
        assert c['original_input']==prior.mesh_digest(mesh)
        assert c['parent']==(-1 if i==0 else 0) and c['depth']==(0 if i==0 else 1)
        assert c['reversed'] is (sign<0)
        if sign<0:mesh['faces']=[[f[0],f[2],f[1]] for f in mesh['faces']]
        assert c['checker_input']==prior.mesh_digest(mesh)
        prior.check_region(c['region'],mesh,1,0)
        meshes.append(mesh)
    assert row['face_components']==labels
    assert sorted(j for c in row['components'] for j in c['vertices'])==list(range(len(material['vertices'])))
    combinations=list(itertools.combinations(range(count),2))
    assert row['pair_count']==len(row['pairs'])==len(combinations)
    distances=[];passed=[]
    for pair,(a,b) in zip(row['pairs'],combinations):
        assert pair['components']==[a,b] and pair['completed'] is True
        assert pair['applicable'] is True and pair['contact'] is False
        ca,cb=row['components'][a],row['components'][b]
        ma,mb=meshes[a],meshes[b]
        assert pair['a_input']==ca['checker_input'] and pair['b_input']==cb['checker_input']
        assert pair['a']==ca['region']['digest'] and pair['b']==cb['region']['digest']
        assert F(pair['minimum_exact'])==F(row['minimum'])
        points=[geo.vec(v) for v in ma['vertices']+mb['vertices']]
        fb=[[j+len(ma['vertices']) for j in f] for f in mb['faces']]
        d=geo.boundary_distance(points,ma['faces'],fb)
        assert d==F(pair['distance']) and d>0
        fa,fb=pair['local_witness']
        assert pair['parent_witness']==[ca['faces'][fa],cb['faces'][fb]]
        assert geo.triangle_distance([geo.vec(ma['vertices'][j]) for j in ma['faces'][fa]],
                                     [geo.vec(mb['vertices'][j]) for j in mb['faces'][fb]])==d
        meets=d>=F(row['minimum'])**2
        assert pair['at_least'] is meets
        distances.append(d);passed.append(meets)
    if not combinations:
        assert row['state']=='NOT_APPLICABLE' and row['applicable'] is False
        assert all(row[k] is None for k in ('contact','at_least','distance','minimum_pair','minimum_components','minimum_local_witness','minimum_parent_witness'))
        accepted=True
    else:
        minimum=min(distances);i=distances.index(minimum);p=row['pairs'][i]
        assert row['applicable'] is True and row['contact'] is False
        assert F(row['distance'])==minimum and row['minimum_pair']==i
        assert row['minimum_components']==p['components']
        assert row['minimum_local_witness']==p['local_witness'] and row['minimum_parent_witness']==p['parent_witness']
        accepted=all(passed);assert row['at_least'] is accepted
    assert row['outcome']==(0 if accepted else 2)
    assert row['published'] is accepted and row['same_input'] is accepted and row['same_binary'] is accepted
    assert row['nef_exceeded'] is False and 0<row['nef_peak']<=row['nef_limit']
    return error


def check_file(path):
    rows=[json.loads(line) for line in Path(path).read_text().splitlines()]
    assert len(rows)==len(IDS) and {r['id'] for r in rows}==IDS
    return rows,max(verify(row) for row in rows)


if __name__=='__main__':
    geo.self_test();binary.analytic_self_test()
    _,error=check_file(sys.argv[1])
    print(f'PASS:7 actual component-gate calls; independent component coverage/mesh hashes/BBW/rational moments/all-pair distances/parent witnesses; max moment error {error:.17g}; no full thickness/seal/paper grant')
