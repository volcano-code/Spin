"""Existing formal benchmark with explicit certified integral selection."""
from pathlib import Path
import sys
from run_core_benchmark import main

if __name__=='__main__':
    main(Path(sys.argv[1]),'small',1,'NATIVE_GMP_SURFACE_V1',
         'EXACT_RATIONAL_FIXED_CLASSIFIED_V1','EVALUATION_CHECKPOINT_V1',
         'LAYERED_SAME_CALL_V1','CONSTRUCTOR_PROOF_REUSE_V1','OWNED_INTERIOR_DIFFERENCE_V1',
         integral_profile='CERTIFIED_MPFR256_GMP_FALLBACK_V1')
