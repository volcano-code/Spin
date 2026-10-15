"""Frozen image-derived Flower request; no solver answers or cut selection."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys

RIM_WIDTH = .015
LENGTH_METERS = .020
PLA_DENSITY = 1200.
STEEL_DENSITY = 7800.


def polygons(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'Polygon':
        return [geometry]
    if geometry.geom_type == 'MultiPolygon':
        return list(geometry.geoms)
    raise ValueError('UNSUPPORTED_PLANAR_COMPONENT: ' + geometry.geom_type)


def boundary_loops(geometry):
    from shapely.geometry.polygon import orient
    result = []
    for polygon in sorted(polygons(geometry), key=lambda p: (*p.bounds, p.area)):
        polygon = orient(polygon, sign=1.)
        result.append([list(p) for p in polygon.exterior.coords[:-1]])
        result.extend([list(p) for p in ring.coords[:-1]] for ring in polygon.interiors)
    return result


def circle(radius, n=96):
    from shapely.geometry import Polygon
    return Polygon([(radius*math.cos(2*math.pi*i/n),
                     radius*math.sin(2*math.pi*i/n)) for i in range(n)])


def exact_profile_moments(loops):
    """Exact rational Green integrals of the stored binary64 polygon vertices."""
    total = [F(0) for _ in range(6)]
    for loop in loops:
        for a, b in zip(loop, loop[1:]+loop[:1]):
            x,y,u,v = map(F, (*a,*b))
            c=x*v-u*y
            row=[c/2,c*(x+u)/6,c*(y+v)/6,
                 c*(x*x+x*u+u*u)/12,c*(y*y+y*v+v*v)/12,
                 c*(2*x*y+x*v+u*y+2*u*v)/24]
            total=[s+t for s,t in zip(total,row)]
    return total


def exact_extrusion_moments(layer):
    a,x,y,xx,yy,xy=exact_profile_moments(layer['loops'])
    low,high=map(F,layer['z_interval'])
    rho=F(layer['volume_density']); h=high-low
    z=(high*high-low*low)/2; zz=(high**3-low**3)/3
    return [rho*q for q in [a*h,x*h,y*h,a*z,xx*h,yy*h,a*zz,y*z,x*z,xy*h]]


def build_fixture(source):
    from shapely.geometry import Polygon
    data=json.loads((source/'preview-input.json').read_text())
    original=Polygon(data['loops'][0],data['loops'][1:])
    if not original.is_valid or len(original.interiors)!=48:
        raise ValueError('INVALID_OR_CHANGED_SOURCE_48_HOLES')
    bore=Polygon(data['loops'][data['bearing_loop_index']])
    outer_race_outside=circle(.128)
    if not bore.covers(outer_race_outside):
        raise ValueError('OUTER_RACE_DOES_NOT_FIT_SOURCE_BORE')
    seat=bore.difference(outer_race_outside)
    outer_race=outer_race_outside.difference(circle(.104))
    # True metric erosion protects every original boundary. Scaling would
    # move holes and may cross slender ligaments in this 48-hole geometry.
    design=original.buffer(-RIM_WIDTH, join_style=2)
    fixed_middle=original.difference(design)
    if design.is_empty or not design.is_valid or not original.covers(design):
        raise ValueError('INVALID_ERODED_DESIGN')
    original_loops=data['loops']
    layers=[dict(id='bottom',loops=original_loops,z_interval=[-.08,-.06],volume_density=1.),
            dict(id='top',loops=original_loops,z_interval=[.06,.08],volume_density=1.),
            dict(id='fixed-middle-boundary-rim',loops=boundary_loops(fixed_middle),z_interval=[-.06,.06],volume_density=1.),
            dict(id='bearing-seat',loops=boundary_loops(seat),z_interval=[-.065,.065],volume_density=1.),
            dict(id='outer-race',loops=boundary_loops(outer_race),z_interval=[-.065,.065],volume_density=STEEL_DENSITY/PLA_DENSITY)]
    # Native adapter accepts geometric layers only; IDs stay in our manifest.
    problem=dict(design_region_loops=boundary_loops(design), design_areal_density=.12,
                 extrusion_interval=[-.06,.06],extrusion_volume_density=1.,
                 fixed_extrusion_layers=[{k:v for k,v in layer.items() if k!='id'} for layer in layers],
                 maximum_iterations=100,maximum_function_evaluations=60,
                 residual_tolerance=1e-10,initial_trust_radius=.25,maximum_trust_radius=1.)
    identity=hashlib.sha256((source/'preview-input.json').read_bytes()).hexdigest()
    geometries={'original':original,'design':design,'fixed_middle':fixed_middle,'seat':seat,'outer_race':outer_race}
    admission=dict(source_valid=original.is_valid,original_hole_count=48,
                   design_valid=design.is_valid,design_component_count=len(polygons(design)),
                   design_area=design.area,fixed_middle_area=fixed_middle.area,
                   design_rim_overlap_area=design.intersection(fixed_middle).area,
                   middle_partition_symmetric_difference_area=original.symmetric_difference(design.union(fixed_middle)).area,
                   seat_flower_overlap_area=seat.intersection(original).area,
                   race_seat_overlap_area=outer_race.intersection(seat).area,
                   minimum_design_boundary_clearance=design.boundary.distance(original.boundary),
                   checks_scope='PLANAR_GEOMETRY_AND_POLYGONAL_EXTRUSIONS_NOT_MANUFACTURING_OR_OPTIMIZATION')
    if admission['design_rim_overlap_area']>0 or admission['middle_partition_symmetric_difference_area']>0:
        raise ValueError('INVALID_MIDDLE_PARTITION')
    fixed=[F(0)]*10
    for layer in layers:
        fixed=[a+b for a,b in zip(fixed,exact_extrusion_moments(layer))]
    record=dict(schema='paper-flower-method-input-v1',classification='RECONSTRUCTED_METHOD_VALIDATION',
                source_preview_sha256=identity,source=str(source),rim_width=RIM_WIDTH,
                rim_rule='metric inward erosion of all outer and hole boundaries; mitre joins',
                retained_original_loops=original_loops,bearing_loop_index=data['bearing_loop_index'],
                meters_per_model_unit=LENGTH_METERS,pla_density_kg_m3=PLA_DENSITY,
                outer_race_density_kg_m3=STEEL_DENSITY,layers=layers,
                fixed_moments_order=['M','LX','LY','LZ','QXX','QYY','QZZ','QYZ','QXZ','QXY'],
                fixed_exact_rational_moments=[str(x) for x in fixed],fixed_model_mass_moments=[float(x) for x in fixed],
                excluded_from_rotor=['upper_cap','lower_cap','shaft','inner_race'],
                bearing_rolling_elements='kinematic bearing constraint; rolling-element mass not modeled',
                target_cut_position=None,postsolve_rotation=False,
                protocol=dict(seeds=[0,1],maximum_attempts=8,maximum_solver_evaluations=60),
                admission=admission)
    return problem,record,geometries


def freeze(source,out):
    from solid_view import extrude
    import trimesh
    problem,record,geometries=build_fixture(source)
    # Validate all mesh exports before creating the immutable output directory.
    meshes={}
    layer_geometries=[geometries['original'],geometries['original'],geometries['fixed_middle'],geometries['seat'],geometries['outer_race']]
    for layer,geometry in zip(record['layers'],layer_geometries):
        # Group holes by their actual parent Polygon. A point-in-outer-ring
        # heuristic misassigns island boundaries inside another polygon's hole.
        meshes[layer['id']]=trimesh.util.concatenate([
            extrude(boundary_loops(polygon),*layer['z_interval']) for polygon in polygons(geometry)])
    meshes['original-flower']=extrude(record['retained_original_loops'],-.08,.08)
    meshes['design']=extrude(problem['design_region_loops'],-.06,.06)
    out.mkdir(parents=True,exist_ok=False)
    meshchecks={}
    for name,mesh in meshes.items():
        mesh.apply_scale(LENGTH_METERS)
        path=out/(name+'.obj');mesh.export(path,digits=17)
        loaded=trimesh.load(path,force='mesh',process=False)
        if not loaded.is_watertight or not loaded.is_winding_consistent or loaded.volume<=0:
            raise ValueError('OBJ_RELOAD_INVALID: '+name)
        error=abs(loaded.volume-mesh.volume)/mesh.volume
        if error>1e-10:raise ValueError('OBJ_VOLUME_CHANGED: '+name)
        meshchecks[name]=dict(watertight=True,winding_consistent=True,positive_volume=True,
                              reload_relative_volume_error=error,volume_m3=loaded.volume,
                              sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    record['obj_checks']=meshchecks
    manifest=dict(schema='faster-planar-method-fixture-v1',scope='RECONSTRUCTED_METHOD_VALIDATION',
                  cases=[dict(id='paper-image-derived-flower-v1',mode='FIDGET_SPINNER',asset_class='RECONSTRUCTED',planar_problem=problem)])
    for name,value in [('flower-cases.json',manifest),('input.json',record)]:
        (out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    (out/'input-sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file()},indent=2)+'\n')
    print(json.dumps(record['admission'],indent=2))


if __name__=='__main__':
    freeze(Path(sys.argv[1]),Path(sys.argv[2]))
