"""Compile only observer propagation in an isolated native overlay."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]


def main(output,reuse=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    commands=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())
    compile_template=commands[0]
    link=list(commands[1])
    if reuse is not None:
        receipt=json.loads((Path(reuse)/'receipt.json').read_text())
        link=list(receipt['commands'][-1]['argv'])
    receipts=[]
    started=time.monotonic()
    with (output/'build.log').open('x') as log:
        stems=('QuadricMomentEvaluator',) if reuse is not None else (
            'PolyhedralMomentDifferential','QuadricMomentEvaluator','FullRankMomentEvaluator','CanonicalPolyhedralMomentDifferential')
        for stem in stems:
            source_stem='PolyhedralMomentDifferential' if stem=='CanonicalPolyhedralMomentDifferential' else stem
            source=ROOT/f'src/spin/{source_stem}.cpp'
            obj=output/f'{stem}.cpp.o'
            argv=list(compile_template)
            if stem=='FullRankMomentEvaluator':
                originals=json.loads((ROOT/'build/heart-exact-input-20261011-v1/commands.json').read_text())
                argv=list(next(c for c in originals if '-c' in c and c[c.index('-c')+1].endswith('/FullRankMomentEvaluator.cpp')))
            if stem=='CanonicalPolyhedralMomentDifferential':
                argv.extend(('-DdifferentiatePolyhedralMaterialMoments=differentiateCanonicalPolyhedralMaterialMoments',
                    '-DdifferentiateRestrictedPolyhedralMaterialMoments=differentiateCanonicalRestrictedPolyhedralMaterialMoments',
                    '-DdifferentiateRationalPolyhedralMaterialMoments=differentiateCanonicalRationalPolyhedralMaterialMoments'))
            for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
                argv[argv.index(flag)+1]=value
            begin=time.monotonic()
            subprocess.run(argv,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=1800)
            receipts.append(dict(argv=argv,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest(),seconds=time.monotonic()-begin))
            for index,entry in enumerate(link):
                if entry.endswith('/'+stem+'.cpp.o'):link[index]=str(obj)
            if stem=='CanonicalPolyhedralMomentDifferential':
                if not any(x.endswith('/'+stem+'.cpp.o') for x in commands[1]) and reuse is None:
                    link.insert(next(i for i,x in enumerate(link) if x.endswith('.a')),str(obj))
        target=output/'_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        link[link.index('-o')+1]=str(target)
        subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=300)
        receipts.append(dict(argv=link,retained_input_hashes={entry:hashlib.sha256(Path(entry).read_bytes()).hexdigest()
            for entry in link if entry.endswith(('.o','.a')) and Path(entry).is_file()}))
    with (output/'receipt.json').open('x') as stream:
        json.dump(dict(schema='native-overlay-build-v1',commands=receipts,seconds=time.monotonic()-started,
            native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),production_admission=False),stream)


if __name__=='__main__': main(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else None)
