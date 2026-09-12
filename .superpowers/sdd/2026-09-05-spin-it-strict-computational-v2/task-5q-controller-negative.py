"""Check rejection of corrupted actual source-moment exports."""
import pathlib
import subprocess
import sys
import tempfile

checker=pathlib.Path(__file__).with_name('task-5q-controller-check.py')
original=pathlib.Path(sys.argv[1]).read_text().splitlines()
starts=[i for i,line in enumerate(original) if line.startswith('CASE ')]
assert len(starts)==2
with tempfile.TemporaryDirectory(prefix='task5q-negative-') as directory:
    target=pathlib.Path(directory)/'corrupt.txt'
    count=0
    for case in range(2):
        begin=starts[case]
        end=next(i for i in range(begin,len(original)) if original[i]=='END')
        for prefix in ('JB ','JH ','MOVED ','MOMENTS ','VALUEMOMENTS ','FLAGS '):
            lines=original.copy()
            i=next(i for i in range(begin,end) if lines[i].startswith(prefix))
            row=lines[i].split()
            field=2 if prefix in ('JB ','JH ','MOVED ') else 1
            row[field]='0' if prefix=='FLAGS ' else str(float(row[field])+.01)
            lines[i]=' '.join(row)
            target.write_text('\n'.join(lines)+'\n')
            result=subprocess.run([sys.executable,str(checker),str(target)],capture_output=True,text=True)
            assert result.returncode!=0 and 'AssertionError' in result.stderr,(case,prefix,result.stderr)
            count+=1
    print(f'PASS: {count} deliberate derivative/geometry/moment/mode corruptions rejected.')
