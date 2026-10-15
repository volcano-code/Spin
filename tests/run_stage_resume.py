"""Repeatable isolated directory for the native in-process stage regression."""
from pathlib import Path
import subprocess
import sys
import tempfile

if len(sys.argv)!=2:raise ValueError('native stage-resume executable required')
with tempfile.TemporaryDirectory(prefix='spin-stage-resume-') as temporary:
    subprocess.run([str(Path(sys.argv[1]).resolve(strict=True)),
                    str(Path(temporary)/'attempts')],check=True,timeout=180)
