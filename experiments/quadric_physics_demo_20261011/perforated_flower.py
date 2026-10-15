"""Pre-registered reconstructed Fig.9-style structure, not author geometry.

Only builds inputs. No solver output, initial answer, or selected cut is read.
All lengths are model units; SI conversion is recorded separately.
"""
from pathlib import Path
import hashlib
import json
import math
import sys

def build_input():
    n = 240
    angles = [2*math.pi*i/n for i in range(n)]
    outer = [[(1+.12*math.cos(12*t))*math.cos(t),
              (1+.12*math.cos(12*t))*math.sin(t)] for t in angles]
    inner = [[.94*x, .94*y] for x,y in outer]
    def circle(radius):
        return [[radius*math.cos(2*math.pi*i/64), radius*math.sin(2*math.pi*i/64)]
                for i in range(64)]
    bore, bearing = circle(.12), circle(.22)
    holes = []
    for k in range(12):
        angle = k*2*math.pi/12
        # Fixed smooth asymmetry creates a genuine COM correction problem.
        radial = .255*(1+.16*math.cos(angle-.4))
        tangential = .068*(1+.12*math.sin(angle+.2))
        loop = []
        for j in range(32):
            t = 2*math.pi*j/32
            x, y = .65+radial*math.cos(t), tangential*math.sin(t)
            loop.append([x*math.cos(angle)-y*math.sin(angle),
                         x*math.sin(angle)+y*math.cos(angle)])
        holes.append(loop[::-1])
    original = [outer, bore[::-1], *holes]
    problem = dict(
        design_region_loops=[inner, bearing[::-1], *holes],
        design_areal_density=.14, extrusion_interval=[-.07,.07],
        extrusion_volume_density=1.,
        fixed_extrusion_layers=[
            dict(loops=original, z_interval=[-.10,-.07], volume_density=1.),
            dict(loops=original, z_interval=[.07,.10], volume_density=1.),
            dict(loops=[outer,inner[::-1]], z_interval=[-.07,.07], volume_density=1.),
            dict(loops=[bearing,bore[::-1]], z_interval=[-.07,.07], volume_density=1.)],
        maximum_iterations=100, maximum_function_evaluations=60,
        residual_tolerance=1e-10, initial_trust_radius=.25, maximum_trust_radius=1.)
    geometry = dict(original_loops=original, outer=outer, inner=inner,
        holes=holes, bearing=bearing, bore=bore, petal_count=12,
        meters_per_model_unit=.020, density_kg_m3=1200,
        classification='RECONSTRUCTED_METHOD_VALIDATION',
        author_original=False, geometric_parameters_are_author_values=False,
        reference='https://pub.ista.ac.at/group_wojtan/projects/2024_Hafner_SIF/sif-final.pdf',
        reference_locator='Figure 9 / Section 6.3; structural reference only',
        structure='fixed perforated top and bottom; optimized middle; fixed rim and bearing ring')
    return problem, geometry

def freeze(directory):
    problem, geometry = build_input()
    directory.mkdir(parents=True, exist_ok=False)
    manifest = dict(schema='faster-planar-method-fixture-v1',
        scope='RECONSTRUCTED_METHOD_VALIDATION',
        protocol='seeds 0,1; 8 attempts; 60 solver evaluations each; input-only initialization',
        cases=[dict(id='reconstructed-perforated-three-layer-flower-v1',
                    mode='FIDGET_SPINNER', asset_class='RECONSTRUCTED', planar_problem=problem)])
    for name, value in [('flower-cases.json',manifest), ('flower-geometry.json',geometry)]:
        (directory/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    identities = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in directory.iterdir() if p.is_file()}
    (directory/'input-sha256.json').write_text(json.dumps(identities,indent=2)+'\n')
    print(directory)

if __name__ == '__main__':
    freeze(Path(sys.argv[1]))
