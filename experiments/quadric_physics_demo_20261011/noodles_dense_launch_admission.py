"""Explicit dense-v2 wrapper; leaves the v1 launch CLI and semantics available."""
from noodles_launch_admission import main

if __name__=='__main__':raise SystemExit(main(input_version='dense-v2'))
