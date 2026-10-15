"""Persist exact zero-area witnesses from a failed derived OBJ batch."""
from pathlib import Path
import hashlib
import json
import os
import sys
import trimesh
from paper_flower_output import exact_zero_area_faces


def diagnose(folder,report):
    record=json.loads((folder/'output-preparation.json').read_text())
    rows=[]
    for name,info in record['files'].items():
        path=folder/info['file']
        mesh=trimesh.load(path,force='mesh',process=False)
        witnesses=exact_zero_area_faces(mesh)
        if witnesses:rows.append(dict(material_or_source_role=name,source_file=str(path),
                                      source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                      exact_zero_area_face_count=len(witnesses),witnesses=witnesses))
    result=dict(schema='paper-flower-obj-exact-zero-area-witness-v1',
                source_preparation_sha256=hashlib.sha256((folder/'output-preparation.json').read_bytes()).hexdigest(),
                records=rows,definition='Exact rational cross product of stored binary64 triangle coordinates is zero.',
                remedy='Separate derived representation: binary64 round-trip OBJ serialization and constrained planar cap retriangulation. Frozen solver input and solved circle unchanged.')
    with report.open('x') as stream:
        stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
    print(json.dumps({row['material_or_source_role']:row['exact_zero_area_face_count'] for row in rows}))


if __name__=='__main__':diagnose(Path(sys.argv[1]),Path(sys.argv[2]))
