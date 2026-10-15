"""Run the exact captured Boolean operands, with no optimizer or fallback."""
import hashlib
import json
from pathlib import Path
import sys
import time


def main(path, output):
    import _cherchi_core_boolean as bridge
    body = json.loads(path.read_text())
    sha = body.pop('operand_packet_sha256')
    if hashlib.sha256(json.dumps(body, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest() != sha:
        raise ValueError('CAPTURE_PACKET_CHANGED')
    expected = body['context']['backend_identity']
    actual = 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:' + hashlib.sha256(Path(bridge.__file__).read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError('CAPTURE_BACKEND_CHANGED')
    print('CAPTURE_VERIFIED', sha, flush=True)
    start = time.monotonic()
    result = bridge.boolean(*(body[name] for name in
        ('design_vertices','design_faces','source_vertices','source_faces','operation')))
    print('BOOLEAN_RETURNED', time.monotonic()-start, flush=True)
    # Preserve returned exact output for a later independent material audit.
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        json.dump(result, stream, sort_keys=True, allow_nan=False)


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
