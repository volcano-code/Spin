"""Validate predecessor bindings without reusing job-reserved key names."""
from pathlib import Path
import hashlib


def inherited_files(identities):
    files = {}
    for key, identity in identities.items():
        path = Path(identity['path'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != identity['sha256']:
            raise ValueError('PREPARATION_DEPENDENCY_CHANGED:' + str(path))
        files['preparation_dependency:' + key] = path
    return files
