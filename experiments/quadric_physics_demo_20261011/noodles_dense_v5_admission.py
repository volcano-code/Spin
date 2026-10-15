"""Exact v5 source/union admission; no solver or simulation invocation."""
import argparse
from pathlib import Path
from noodles_admission import admit

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    admit(args.input,args.output,input_version='dense-v5')
