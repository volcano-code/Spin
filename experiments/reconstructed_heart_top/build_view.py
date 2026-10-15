"""Build the existing strict viewer against the experiment's actual extension.

Reuses installed OpenGL libraries, never replaces the old executable or module.
"""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

root=Path(__file__).resolve().parents[2]
old=root/'build/task-7v1-viewer'
out=root/'build/reconstructed-heart-view-20260913'
out.mkdir(exist_ok=False)
flags={}
for line in (old/'CMakeFiles/make_it_stand_viewer.dir/flags.make').read_text().splitlines():
    if line.startswith('CXX_'):
        k,v=line.split(' = ',1); flags[k]=shlex.split(v.replace('\\"','"'))
defines=[v for v in flags['CXX_DEFINES'] if not v.startswith('-DSEALED_SPIN_')]
defines += ['-DSEALED_SPIN_PYTHON="/home/lin/anaconda3/bin/python3.13"',
    f'-DSEALED_SPIN_SOURCE_PYTHON="{root}/python"',
    f'-DSEALED_SPIN_EXTENSION="{root}/build/task-5c-no-faster-20260909/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so"']
sources=[root/'src/viewer/SealedSpinView.cpp',Path(__file__).with_name('sealed_view_main.cpp')]
commands=[]; objects=[]
for i,src in enumerate(sources):
    obj=out/f'view{i}.o'; objects.append(str(obj))
    command=['/usr/bin/c++',*defines,*flags['CXX_INCLUDES'],'-O1','-DNDEBUG','-std=c++17',
             '-fPIC','-c',str(src),'-o',str(obj)]
    commands.append(command); subprocess.run(command,check=True,cwd=root)
libs=[old/'libimguizmo.a',old/'imgui/libimgui.a',old/'glad/libglad.a',old/'glfw/src/libglfw3.a',
      Path('/home/lin/anaconda3/lib/libpython3.13.so'),
      Path('/usr/lib/x86_64-linux-gnu/libGLX.so'),Path('/usr/lib/x86_64-linux-gnu/libOpenGL.so')]
command=['/usr/bin/c++',*objects,*map(str,libs),'-Wl,-rpath,/usr/lib/x86_64-linux-gnu:/home/lin/anaconda3/lib',
         '-lrt','-lm','-ldl','-pthread','-o',str(out/'sealed-heart-view')]
commands.append(command); subprocess.run(command,check=True,cwd=root)
record=dict(commands=commands,files={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [*sources,*libs,out/'sealed-heart-view']},
            scope='Existing strict read-only OpenGL viewer; no optimization or physics in viewer')
(out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
print(out/'sealed-heart-view')
