#!/usr/bin/env python

import shutil, os, argparse, sys, stat
sys.path.append("cppSetupTools/scripts/pyUtils")
sys.path.append("cppSetupTools/scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper
def main():
    name = "bhtsne"
    libs = "armadillo:8.200.0,libpca:develop,bibcpp:develop"
    args = genHelper.parseNjhConfigureArgs()
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    Utils.run(cmd)
    
main()


