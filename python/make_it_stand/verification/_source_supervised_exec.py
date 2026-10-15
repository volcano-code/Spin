"""Internal POSIX exec shim: set per-process address-space cap before exec."""
import os
import resource
import sys


if __name__=='__main__':
    memory=int(sys.argv[1])
    resource.setrlimit(resource.RLIMIT_AS,(memory,memory))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    os.execv(sys.argv[2],sys.argv[2:])
