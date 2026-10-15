"""One-time exact represented-entity comparison, outside normal timing."""
import argparse
import json
from pathlib import Path
import numpy as np


def main(inputs,output):
    import _make_it_stand as native
    from make_it_stand.verification.faster_material import read_exact_boolean_sidecar_coordinates
    from make_it_stand.verification.spin_problem import MeshSnapshot
    from make_it_stand.verification.cherchi_motion_differential import normalize_closed_component_orientation
    import sys
    root=Path(__file__).resolve().parents[2]
    sys.path.insert(0,str(root/'experiments/spin_it_reproduction/strict'))
    from run_faster_reference_boolean import _read_obj
    output.mkdir(parents=True,exist_ok=False)
    records=[]
    for state in ('ellipsoid','one_sheet','two_sheet'):
        folder=inputs/state
        dv,df=_read_obj(folder/'design.obj');sv,sf=_read_obj(folder/'source.obj')
        evaluated=json.loads((inputs/(state+'-cgal-evaluation.json')).read_text())
        operation='intersection' if evaluated['quadric_solid_is_positive_region'] else 'difference'
        row=dict(state=state,operation=operation,optimization_admission=False)
        try:
            construction=native.libigl_exact_boolean_construction(dv,df,sv,sf,operation)
            if not construction['succeeded']:raise ValueError('CGAL_EXACT_CONSTRUCTION_FAILED')
            av,af=_read_obj(folder/'cherchi-0.obj')
            rational,_=read_exact_boolean_sidecar_coordinates(av,Path(str(folder/'cherchi-0.obj')+'.exact.tsv'))
            snapshot=MeshSnapshot.from_rationals(rational,af)
            cv,cf=snapshot.arrays()
            cf,_=normalize_closed_component_orientation(cv,cf,rational_vertices=snapshot.rational_vertices)
            a=native.make_exact_material_region_from_rationals(construction['vertices'],construction['faces'])
            b=native.make_exact_material_region_from_rationals(rational,cf)
            row['cgal_region_valid']=a.facts['valid'];row['cherchi_region_valid']=b.facts['valid']
            row['relation']=dict(native.relate_exact_material_regions(a,b))
        except Exception as exc:
            row.update(exception_type=type(exc).__name__,exception=str(exc))
        records.append(row)
        with (output/(state+'.json')).open('x') as stream:json.dump(row,stream)
    with (output/'report.json').open('x') as stream:json.dump(dict(records=records),stream)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();main(args.inputs,args.output)
