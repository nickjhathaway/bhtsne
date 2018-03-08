#!/usr/bin/env python3



import subprocess, sys, os, argparse,shutil
from collections import namedtuple, defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "cppSetupTools"))
from setup import *


if __name__ == '__main__':
  SetupRunner.runSetup()
    
    
    
