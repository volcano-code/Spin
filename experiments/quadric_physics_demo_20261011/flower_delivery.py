"""Publish final-only assembly separately from diagnostic material volumes."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def final_parts(record):
    names=['retained',*record['rotor_layers'],*record['clamp_layers']]
    forbidden={'design','cavity','original-flower','original-middle','bearing-bore-region'}
    if forbidden.intersection(names):
        raise ValueError('DIAGNOSTIC_SOLID_IN_FINAL_ASSEMBLY')
    if len(set(names))!=len(names) or any(name not in record['files'] for name in names):
        raise ValueError('FINAL_ASSEMBLY_ROLE_OR_FILE_MISMATCH')
    return names


def publish(folder, material, saved_check, output):
    from flower_run_physics import admit
    folder,material,saved_check,output=map(Path,(folder,material,saved_check,output))
    record,proof=admit(material,saved_check.parent,folder,saved_material_check=saved_check)
    names=final_parts(record)
    assembly=output/'final-assembly'
    diagnostic=output/'diagnostics'
    output.mkdir(parents=True,exist_ok=False)
    assembly.mkdir();diagnostic.mkdir()
    copied={}
    for name,info in record['files'].items():
        source=material/info['file']
        if sha(source)!=info['sha256']:
            raise ValueError('MATERIAL_FILE_IDENTITY_MISMATCH: '+name)
        destination=(assembly if name in names else diagnostic)/info['file']
        shutil.copyfile(source,destination)
        copied[str(destination.relative_to(output))]=sha(destination)
    properties=dict(schema='flower-final-assembly-physical-properties-v1',
        units='meters',classification='RECONSTRUCTED_METHOD_VALIDATION',
        original=record['original'],optimized=record['optimized'],clamp=record['clamp'],
        rotor_parts=['retained',*record['rotor_layers']],fixed_clamp_parts=list(record['clamp_layers']),
        pla_density_kg_m3=1200.,
        material_density_ratios={'retained':1.,**record['rotor_layers'],**record['clamp_layers']},
        independent_balance=proof['independent_saved_balance'],manufacturing_qualified=False)
    (assembly/'physical-properties.json').write_text(json.dumps(properties,indent=2)+'\n')
    (assembly/'README.txt').write_text(
        'Final assembly only: load every OBJ in this directory together.\n'
        'retained.obj is the optimized middle layer. top.obj and bottom.obj are fixed covers.\n'
        'For an exploded view, move the covers only in display space.\n'
        'Units: meters. Diagnostic cavity/design/original volumes are in a different package.\n'
        'Clamps are fixed components; do not include them in rotor inertia.\n')
    for subdir,archive in [(assembly,'final-assembly-obj.zip'),(diagnostic,'diagnostic-obj.zip')]:
        with zipfile.ZipFile(output/archive,'x',compression=zipfile.ZIP_DEFLATED) as stream:
            for path in sorted(subdir.iterdir()):
                stream.write(path,str(path.relative_to(output)))
    manifest=dict(schema='flower-separated-delivery-v1',
        material_preparation_sha256=sha(material/'output-preparation.json'),
        saved_check_sha256=sha(saved_check),
        final_parts=names,files=copied,
        archives={name:sha(output/name) for name in ['final-assembly-obj.zip','diagnostic-obj.zip']})
    (output/'delivery-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for key in ('folder','material','saved-check','output'):
        parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args()
    publish(args.folder,args.material,args.saved_check,args.output)
