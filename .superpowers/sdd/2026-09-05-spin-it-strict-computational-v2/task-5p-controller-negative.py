"""Reject deliberately corrupted finite-fixture exports; no project imports."""
import pathlib
import subprocess
import sys
import tempfile

checker = pathlib.Path(__file__).with_name('task-5p-controller-check.py')
original = pathlib.Path(sys.argv[1]).read_text().splitlines()

def mutate(case, kind):
    lines = original.copy()
    starts = [i for i, line in enumerate(lines) if line.startswith('CASE ')]
    begin = starts[case]
    end = next(i for i in range(begin, len(lines)) if lines[i] == 'END')
    prefix = {'missing_ref': 'REF ', 'triangle_flip': 'TRI ',
              'wrong_owner': 'PATCH ', 'wrong_moment': 'MOMMOVED ',
              'wrong_motion': 'MOVED ', 'duplicate_key': 'KEY '}[kind]
    i = next(i for i in range(begin, end) if lines[i].startswith(prefix))
    row = lines[i].split()
    if kind == 'missing_ref':
        del lines[i]
        return lines
    if kind == 'triangle_flip': row[3], row[4] = row[4], row[3]
    if kind == 'wrong_owner': row[7] = row[8]
    if kind in ('wrong_moment', 'wrong_motion'): row[2] = str(float(row[2]) + .01)
    if kind == 'duplicate_key':
        following = next(j for j in range(i+1, end) if lines[j].startswith(prefix))
        row[2:] = lines[following].split()[2:]
    lines[i] = ' '.join(row)
    return lines

with tempfile.TemporaryDirectory(prefix='task5p-negative-') as directory:
    target = pathlib.Path(directory) / 'corrupt.txt'
    count = 0
    for case in range(3):
        for kind in ('missing_ref', 'triangle_flip', 'wrong_owner', 'wrong_moment', 'wrong_motion', 'duplicate_key'):
            target.write_text('\n'.join(mutate(case, kind)) + '\n')
            result = subprocess.run([sys.executable, str(checker), str(target)], capture_output=True, text=True)
            assert result.returncode != 0 and 'AssertionError' in result.stderr, (case, kind, result.stdout, result.stderr)
            count += 1
    print(f'PASS: {count} corruptions rejected across 3 fixtures (refs, winding, ownership, moments, motion, keys).')
