#!/usr/bin/env python



import subprocess, sys, os, argparse,shutil
from collections import namedtuple
sys.path.append(os.path.join(os.path.dirname(__file__), "scripts/pyUtils"))
sys.path.append(os.path.join(os.path.dirname(__file__), "scripts/setUpScripts"))
from utils import Utils
from genFuncs import genHelper 
from color_text import ColorText as CT
import pickle, datetime

#tuples
BuildPaths = namedtuple("BuildPaths", 'url build_dir build_sub_dir local_dir')
LibNameVer = namedtuple("LibNameVer", 'name version')
GitRefs = namedtuple("GitRefs", "branches tags")


class LibDirMaster():
    def __init__(self,externalLoc):
        self.base_dir = os.path.abspath(externalLoc); #top dir to hold tars,build, local directories
        
        self.ext_tars = os.path.join(self.base_dir, "tarballs") #location to keep tarballs of programs/libraries downloads
        self.ext_build = os.path.join(self.base_dir, "build") #location for the building of programs/libraries
        self.install_dir = os.path.join(self.base_dir, "local") #location for the final install of programs/libraries
        self.cache_dir = os.path.join(self.base_dir, ".cache")
        
        Utils.mkdir(self.ext_tars) #tar storage directory
        Utils.mkdir(self.ext_build) #build directory
        Utils.mkdir(self.install_dir) #local directory
        Utils.mkdir(self.cache_dir) #cache directory

def joinNameVer(libNameVerTup):
    return os.path.join(libNameVerTup.name, libNameVerTup.version, libNameVerTup.name)




class CPPLibPackageVersion():
    def __init__(self, name, version, bPaths, depends):
        self.nameVer_ = LibNameVer(name, version)
        self.depends_ = depends #should be a list of LibNameVer
        self.bPaths_ = bPaths
        self.includePath_ = os.path.join(joinNameVer(self.nameVer_), "include")
        self.additionalIncludeFlags_ = []
        self.additionalIncludePaths_ = []
        self.libPath_ = os.path.join(joinNameVer(self.nameVer_), "lib")
        self.additionalLdFlags_ = []
        self.libName_ = name
        self.altLibName_ = ""
        
        
    def getDownloadUrl(self):
        ret = self.bPaths_.url
        if str(self.bPaths_.url).endswith(".git"):
            ret = self.bPaths_.url.replace(".git","/archive/" + str(self.nameVer_.version) + ".tar.gz").replace("git@github.com:", "https://github.com/")
        return ret
    
    def getIncludeFlags(self, localPath):
        ret = ""
        if(len(self.includePath_) > 0):
            ret = "-isystem" + str(os.path.join(localPath, self.includePath_))
        if len(self.additionalIncludePaths_) > 0:
            for addPath in self.additionalIncludePaths_:
                if len(ret) > 0:
                    ret = ret + " "
                ret = ret + "-isystem" + str(os.path.join(localPath, addPath))
        if len(self.additionalIncludeFlags_) > 0:
            if len(ret)> 0:
                ret = ret + " "
            ret = ret + " ".join(self.additionalIncludeFlags_) 
        return ret
    
    def getLdFlags(self, localPath):
        ret = ""
        retList = []
        libPath = str(os.path.join(localPath,self.libPath_))
        if(len(self.libPath_) > 0):
            retList.append("-Wl,-rpath," + str(libPath))
            retList.append("-L" + str(libPath))
            if len(self.altLibName_) > 0:
                retList.append("-l" + self.altLibName_)
            elif "" != self.libName_:
                retList.append("-l" + self.libName_)
        if len(self.additionalLdFlags_) > 0:
            retList.extend(self.additionalLdFlags_)
        if len(retList) > 0:
            ret = " ".join(retList)                 
        return ret
    

class CPPLibPackage():
    def __init__(self, name, defaultBuildCmd, dirMaster, libType, defaultVersion):
        self.name_ = name
        
        self.defaultVersion_ = defaultVersion.replace("/", "__")
        self.defaultBuildCmd_ = defaultBuildCmd
        self.versions_ = {}
        self.externalLibDir_ = dirMaster
        if "git" != libType and "file" != libType and "git-headeronly" != libType:
            raise Exception("libType should be 'git', 'git-headeronly', or 'file', not " + str(libType))
        self.libType_ = libType #should be git, git-headeronly, or file
        self.bibProject_ = False
    
    def addVersion(self, url, verName, depends=[]):
        verName = verName.replace("/", "__")
        build_dir = os.path.join(self.externalLibDir_.ext_build, self.name_, verName)
        fn = os.path.basename(url)
        fn_noex = fn.replace(".tar.gz", "").replace(".tar.bz2", "").replace(".git", "")
        build_sub_dir = os.path.join(self.externalLibDir_.ext_build, self.name_, verName, self.name_)
        local_dir = os.path.join(self.externalLibDir_.install_dir, self.name_, verName, self.name_)
        self.versions_[verName] = CPPLibPackageVersion(self.name_, verName,BuildPaths(url, build_dir, build_sub_dir, local_dir), depends)
    
    def addHeaderOnlyVersion(self, url, verName, depends=[]):
        '''set up for header only libraries, these just need
         the header copied no need for build_dir build_sub_dir '''
        verName = verName.replace("/", "__")
        local_dir = os.path.join(self.externalLibDir_.install_dir, self.name_, verName, self.name_)
        self.versions_[verName] = CPPLibPackageVersion(self.name_, verName,BuildPaths(url, "", "", local_dir), depends)
        self.versions_[verName].includePath_ = os.path.join(self.name_, verName)
        #self.versions_[verName].includePath_ = joinNameVer(self.versions_[verName].nameVer_)
        self.versions_[verName].libPath_ = ""
        
    def hasVersion(self, version):
        return version in self.versions_
    
    def getVersions(self):
        return sorted(self.versions_.keys())
    
    def getLocalDir(self, version):
        if self.hasVersion(version):
            return self.versions_[version].bPaths_.local_dir
        raise Exception("Error in getLocalDir" + self.name_ + " doesn't have version " + str(version))
    
    def getBuildSubDir(self, version):
        if self.hasVersion(version):
            return self.versions_[version].bPaths_.build_sub_dir
        raise Exception("Error in getBuildSubDir" + self.name_ + " doesn't have version " + str(version))
    
    def getBuildDir(self, version):
        if self.hasVersion(version):
            return self.versions_[version].bPaths_.build_dir
        raise Exception("Error in getBuildDir" + self.name_ + " doesn't have version " + str(version))
    
    def getGitRefs(self, url):
        if not self.libType_.startswith("git"):
            raise Exception("Library " + self.name_ + " is not a git library, type is : " + self.libType_)
        try:
            cap = Utils.runAndCapture("git ls-remote {url}".format(url = url))
        except Exception as inst: 
            try:
                #if the first attempt fail, try doing https instead if that was reason
                url = url.replace("git@github.com:", "https://github.com/")
                cap = Utils.runAndCapture("git ls-remote {url}".format(url = url))
            except Exception as instFallback:
                raise instFallback 
        branches = []
        tags = []
        for line in cap.split("\n"):
            if "" != line:
                lineSplit = line.split()
                if 2 == len(lineSplit):
                    if "heads" in lineSplit[1]:
                        branches.append(lineSplit[1][(lineSplit[1].find("heads/") + 6):])
                    elif "tags" in lineSplit[1] and not lineSplit[1].endswith("^{}"):
                        tags.append(lineSplit[1][(lineSplit[1].find("tags/") + 5):])
        gRefs = GitRefs(branches, tags)
        return (gRefs)
            

class Packages():
    '''class to hold and setup all the necessary paths for 
    downloading, building, and then installing packages/libraries'''
    def __init__(self, externalLoc, args):
        self.dirMaster_ = LibDirMaster(externalLoc); #top dir to hold tars,build, local directories
        self.args = args
        self.packages_ = {} #dictionary to hold path infos
        self.packages_["armadillo"] = self.__armadillo()


    def package(self, name):
        '''get package info if it exists'''
        if name in self.packages_:
            return self.packages_[name]
        raise Exception(name + " not found in paths")

    def __armadillo(self):
        name = "armadillo"
        buildCmd = "mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir} .. && make -j {num_cores} install"
        pack = CPPLibPackage(name, buildCmd, self.dirMaster_, "file", "7.500.2")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-7.500.2.tar.gz", "7.500.2")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-7.400.2.tar.gz", "7.400.2")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-7.300.1.tar.gz", "7.300.1")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-7.100.3.tar.gz", "7.100.3")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-6.700.3.tar.gz", "6.700.3")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-6.200.3.tar.gz", "6.200.3")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-6.100.0.tar.gz", "6.100.0")
        pack.addVersion("http://baileylab.umassmed.edu/sourceCodes/armadillo/armadillo-5.600.2.tar.gz", "5.600.2")
        return pack

    

    
    def getPackagesNames(self):
        return sorted(self.packages_.keys())
    
    def checkForPackVer(self, packVer):
        if packVer.name not in self.packages_:
            raise Exception("Lib " + packVer.name + " not found in libs, options are " + ", ".join(self.getPackagesNames()))
        else:
            if packVer.version.replace("/", "__") not in self.packages_[packVer.name].versions_:
                raise Exception("Version " + packVer.version + " for lib " \
                                + packVer.name + " not found in available versions, options are " \
                                + ", ".join(self.packages_[packVer.name].getVersions()))
        return True
                
    def getLdFlags(self, packVer):
        self.checkForPackVer(packVer)
        return self.packages_[packVer.name].versions_[packVer.version].getLdFlags(self.dirMaster_.install_dir)
    
    def getIncludeFlags(self, packVer):
        self.checkForPackVer(packVer)
        return self.packages_[packVer.name].versions_[packVer.version].getIncludeFlags(self.dirMaster_.install_dir)
    
    def writeMakefile(self, packVers, filename, overwrite = False, append = False):
        if os.path.exists(filename) and not overwrite and not append:
            raise Exception("File: " + str(filename) + " already exists, use --overWrite to overwrite it")
        elif os.path.exists(filename) and overwrite:
            os.remove(filename)
            self.writeMakefile(packVers, filename, overwrite, append)
        elif os.path.exists(filename) and append:
            with open(filename, "a") as f:
                for packVer in packVers:
                    pack = self.package(packVer.name)
                    pvIncFlags = self.getIncludeFlags(packVer)
                    if "" != pvIncFlags:
                        f.write("#" + packVer.name + ":" + packVer.version + " CXXFLAGS\n")
                        f.write("COMLIBS += " + pvIncFlags + "\n")
                    pvLdFlags = self.getLdFlags(packVer)
                    if "" != pvLdFlags:
                        f.write("#" + packVer.name + ":" + packVer.version + " LDFLAGS\n")
                        f.write("LD_FLAGS += " + pvLdFlags + "\n")
                    f.write("\n")
                    f.flush()
        else:
            with open(filename, "a") as f:
                for packVer in packVers:
                    pack = self.package(packVer.name)
                    #if bib project, add the flags of it's dependencies
                    if pack.bibProject_:
                            cmd = "python ./setup.py --compfile compfile.mk --numCores 1 --append --outMakefile {makefileCommon}".format(makefileCommon = os.path.abspath(filename))
                            dir = pack.getBuildSubDir(packVer.version)
                            Utils.run_in_dir(cmd, dir)
                    pvIncFlags = self.getIncludeFlags(packVer)
                    if "" != pvIncFlags:
                        f.write("#" + packVer.name + ":" + packVer.version + " CXXFLAGS\n")
                        f.write("COMLIBS += " + pvIncFlags + "\n")
                    pvLdFlags = self.getLdFlags(packVer)
                    if "" != pvLdFlags:
                        f.write("#" + packVer.name + ":" + packVer.version + " LDFLAGS\n")
                        f.write("LD_FLAGS += " + pvLdFlags + "\n")
                    f.write("\n")
                    f.flush()
    
    def addPackage(self, packVers, packVer):
        packVer = LibNameVer(packVer.name, packVer.version.replace("/", "__"))
        if self.checkForPackVer(packVer):
            pack = self.package(packVer.name)
            for dep in pack.versions_[packVer.version].depends_:
                self.addPackage(packVers, LibNameVer(str(dep.name).lower(), dep.version))
            found = False
            for otherPackVer in packVers:
                if otherPackVer.name == packVer.name:
                    if otherPackVer.version != packVer.version:
                        raise Exception("Version conflict for " + packVer.name + " already have " + otherPackVer.version + " and adding: " + packVer.version)
                    else:
                        found = True
            if not found:
                packVers.append(packVer)
            
                
                
    def isInstalled(self, packVer):
        if os.path.exists(os.path.join(self.dirMaster_.install_dir, joinNameVer(packVer))):
            return True
        else:
            return False
    

    
    def getDefaultIncludeFlags(self):
        return "-I./src/"
    
    def getDefaultLDFlags(self):
        ret = ""
        if Utils.isMac():
            #for dylib path fixing in macs, this gets rid of the name_size limit, which why the hell is there a name size limit
            ret = ret + "-headerpad_max_install_names" 
        return ret


    
    
    
class Setup:
    def __init__(self, args):
        self.extDirLoc = "" # the location where the libraries will be installed
        #if no compile file set up and assume external is next to setup.py
        if not args.compfile:
            self.extDirLoc = "external"
            #self.extDirLoc = os.path.abspath(os.path.join(os.path.dirname(__file__), "external"))
        else:
            self.extDirLoc = os.path.abspath(self.parseForExtPath(args.compfile[0]))
        self.dirMaster_ = LibDirMaster(self.extDirLoc)
        self.args = args # command line arguments parsed by argument parser
        self.setUps = {} # all available set ups
        self.setUpsNeeded = [] # the setups that need to be done
        self.installed = [] # the setups that able to install
        self.failedInstall = [] # the setups that failed
        self.CC = "" # the c compilier being used
        self.CXX = "" # the c++ compilier being used
        self.noInternet_ = False
        if args.noInternet:
            self.noInternet_ = True
        self.__initSetUpFuncs()
        self.__processArgsForCompilers()
        #if we have internet and the cache is more than a day old, clear it
        if Utils.connectedInternet:
            cacheDate = datetime.datetime.fromtimestamp(os.path.getmtime(self.dirMaster_.cache_dir))
            now = datetime.datetime.now()
            if 86400 < (now - cacheDate).total_seconds():
                self.clearCache()
        if args.clearCache:
            self.clearCache()
        self.packages_ = Packages(self.extDirLoc, self.args) # path object to hold the paths for install
        self.__processArgsForSetupsNeeded()
        
    def setup(self):
        if self.args.forceUpdate:
            for set in self.setUpsNeeded:
                if not set.name in self.setUps.keys():
                    print CT.boldBlack( "Unrecognized option ") + CT.boldRed(set.name)
                else:
                    self.rmDirsForLib(set)
                    
        for set in self.setUpsNeeded:
            if not set.name in self.setUps.keys():
                print CT.boldBlack( "Unrecognized option ") + CT.boldRed(set.name)
            else:
                self.__setup(set.name, set.version)

        for p in self.installed:
            print p.name + ":" + str(p.version), CT.boldGreen("installed")

        for p in self.failedInstall:
            print  p.name + ":" + str(p.version), CT.boldRed("failed to install")

    def __initSetUpFuncs(self):
        self.setUps = {
                       "armadillo": self.armadillo

                       }

    def printAvailableSetUps(self):
        self.__initSetUpFuncs()
        print "Available installs:"
        print "To Install use ./setup.py --libs lib1,lib2,lib3"
        print "E.g. ./setup.py --libs armadillo:7.500.2"
        installs = self.setUps.keys()
        installs.sort()
        for set in installs:
            print set
            pack = self.__package(set)
            sys.stdout.write("\t")
            sys.stdout.write(",".join([p.replace("__", "/") for p in pack.getVersions()]))
            sys.stdout.write("\n")
            
    def printGitRefs(self):
        self.__initSetUpFuncs()
        print "Git branches and tags:"
        for set in self.setUpsNeeded:
            print set.name
            pack = self.__package(set.name)
            refs = pack.getGitRefs(pack.versions_[pack.defaultVersion_].bPaths_.url)
            print "\t" + "Branches"
            for b in refs.branches:
                print "\t\t" + b
            print "\t" + "Tags"
            for t in refs.tags:
                print "\t\t" + t

    def __processArgsForSetupsNeeded(self):
        if self.args.libs:
            inLibs = self.args.libs.split(",")
            for lib in inLibs:
                if ":" not in lib.lower():
                    raise Exception("Need to give version for " + lib)
                else:
                    libSplit = lib.split(":")
                    self.packages_.addPackage(self.setUpsNeeded,LibNameVer(libSplit[0].lower(), libSplit[1]))
        if self.args.compfile:
            self.parseSetUpNeeded(self.args.compfile[0])
    
    def __processArgsForCompilers(self):
        if self.args.compfile:
            self.parserForCompilers(self.args.compfile[0])
        # if no compfile need to determine compiler, will default to env CC and CXX
        else:
            self.CC = genHelper.determineCC(self.args)
            self.CXX = genHelper.determineCXX(self.args)
            self.args.CC = self.CC
            self.args.CXX = self.CXX
        if "clang" in self.CXX:
            self.args.clang = True
        else:
            self.args.clang = False

    def parseForExtPath(self, fn):
        args = self.parseCompFile(fn)
        if "EXT_PATH" in args:
            extPath = args["EXT_PATH"].strip()
            extPath = extPath.replace("$(realpath", "")
            extPath = extPath.replace(")", "")
            extPath = extPath.strip()
        else:
            print "did not find external folder location; assuming ./external"
            extPath = "./external"
        return extPath

    def parseSetUpNeeded(self, fn):
        args = self.parseCompFile(fn)
        for k,v in args.iteritems():
            if k.startswith("USE_"):
                if '0' != v:
                    if "#" in v:
                        valSplit = v.split("#")
                        if valSplit[0] == '1':
                            self.packages_.addPackage(self.setUpsNeeded, LibNameVer(k[4:].lower(),valSplit[1]))
                    else:
                        raise Exception("Need to supply version in compfile with USE_PACKAGE#Version")
                

    def parseCompFile(self, fn):
        ret = {}
        with open(fn) as f:
            for line in f:
                if '=' in line:
                    toks = line.split('=')
                    ret[toks[0].strip()] = toks[1].strip()
        return ret

    def parserForCompilers(self, fn):
        args = self.parseCompFile(fn)
        if 'CC' in args:
            self.CC = args['CC']
            self.args.CC = self.CC
        if 'CXX' in args:
            self.CXX = args['CXX']
            self.args.CXX = self.CXX
    
    def rmDirsForLibs(self,libs):
        for l in libs:
            self.rmDirsForLib(l)
    
    def rmDirsForLib(self,packVer):
        if packVer.name not in self.setUps:
            print CT.boldBlack( "Unrecognized package: ") + CT.boldRed(packVer.name)
        else:
            pack = self.__package(packVer.name)
            if not pack.hasVersion(packVer.version):
                raise Exception("No version " + str(packVer.version) + " for " + str(packVer.name))
            p = pack.versions_[packVer.version].bPaths_
            if os.path.exists(p.build_dir):
                print "Removing " + CT.boldBlack(p.build_dir)
                Utils.rm_rf(p.build_dir)
            if os.path.exists(p.local_dir):
                print "Removing " + CT.boldBlack(p.local_dir)
                Utils.rm_rf(p.local_dir)
    

    def __package(self, name):
        return self.packages_.package(name)

    def __setup(self, name, version):
        version = version.replace("/", "__")
        pack = self.__package(name)
        if not pack.hasVersion(version):
            raise Exception("Package " + str(name) + " doesn't have version " + str(version))
        bPath = pack.versions_[version].bPaths_
        if os.path.exists(bPath.local_dir):
            print CT.boldGreen(name + ":" + version), "found at " + CT.boldBlue(bPath.local_dir)
        else:
            print CT.boldGreen(name + ":" + version), CT.boldRed("NOT"), "found; building..."
            try:
                self.setUps[name](version)
                self.installed.append(LibNameVer(name, version))
            except Exception as inst:
                print inst 
                print CT.boldRed("failed to install ") + name + ":" + str(version)
                self.failedInstall.append(LibNameVer(name, version))

    def num_cores(self):
        retCores = Utils.num_cores()
        if self.args.numCores:
            if not self.args.numCores > retCores:
                retCores = self.args.numCores
        else:
            if retCores > 8:
                retCores  = retCores/2
            if 1 != retCores:
                retCores -= 1
        return retCores

    def __buildFromFile(self, packVer, cmd):
        bPath = packVer.bPaths_
        if self.noInternet_:
            newUrl = bPath.url.replace(".git","/archive/" + str(packVer.nameVer_.version) + ".tar.gz").replace("git@github.com:", "https://github.com/")
            bPath = BuildPaths(newUrl, bPath.build_dir, bPath.build_sub_dir, bPath.local_dir)
            base_file = os.path.basename(bPath.url)
            fnp = os.path.join(self.dirMaster_.ext_tars,packVer.nameVer_.name, base_file)
            if not os.path.exists(fnp):
                raise Exception("Could not find file: " + str(fnp))
        else:
            print "\t Getting file..."
            Utils.mkdir(os.path.join(self.dirMaster_.ext_tars, packVer.nameVer_.name))
            fnp = Utils.get_file_if_size_diff(bPath.url, os.path.join(self.dirMaster_.ext_tars, packVer.nameVer_.name))
        Utils.clear_dir(bPath.build_dir)
        Utils.untar(fnp, bPath.build_dir)
        ##probably not the best way to do this as there is no guarantee that there is a directory there
        untaredDir = os.listdir(bPath.build_dir)[0]
        os.rename(os.path.join(bPath.build_dir, untaredDir), bPath.build_sub_dir)
        try:
            Utils.run_in_dir(cmd, bPath.build_sub_dir)
        except:
            print "\t Failed to build, removing {d}".format(d = bPath.local_dir)
            Utils.rm_rf(bPath.local_dir)
            sys.exit(1)
                
    def __buildFromGitBranch(self, packVer, cmd):
        bPath = packVer.bPaths_
        if self.noInternet_:
            self.__buildFromFile(packVer, cmd)
        else:
            if os.path.exists(bPath.build_sub_dir):
                print "pulling from {url}".format(url=bPath.url)
                pCmd = "git checkout " + packVer.nameVer_.version.replace("__", "/") + " && git pull"
                try:
                    Utils.run_in_dir(pCmd, bPath.build_sub_dir)
                except:
                    print "failed to pull from {url} with {cmd}".format(url=bPath.url, cmd = pCmd)
                    sys.exit(1)
            else:
                print "cloning from {url}".format(url=bPath.url)
                cCmd = "git clone -b " + packVer.nameVer_.version.replace("__", "/") + " {url} {d}".format(url=bPath.url, d=bPath.build_sub_dir)
                try:
                    Utils.run(cCmd)
                except:
                    print "failed to clone from {url}".format(url=bPath.url)
                    sys.exit(1)
            try:
                Utils.run_in_dir(cmd, bPath.build_sub_dir)
            except:
                print("Failed to build, removing {d}".format(d = bPath.local_dir))
                Utils.rm_rf(bPath.local_dir)
                sys.exit(1)
    
    def __buildFromGitTag(self, packVer, cmd):
        bPath = packVer.bPaths_
        ##if no internet build from tar file, file needs to be in tarballs folder
        if self.noInternet_:
            self.__buildFromFile(packVer, cmd)
        else:
            if os.path.exists(bPath.build_sub_dir):
                print "pulling from {url}".format(url=bPath.url)
                pCmd = "git checkout master && git pull && git checkout " + packVer.nameVer_.version
                try:
                    Utils.run_in_dir(pCmd, bPath.build_sub_dir)
                except Exception, e:
                    print e
                    print "failed to pull from {url}".format(url=bPath.url)
                    sys.exit(1)
            else:
                print "cloning from {url}".format(url=bPath.url)
                cCmd = "git clone {url} {d}".format(url=bPath.url, d=bPath.build_sub_dir)
                tagCmd = "git checkout {tag}".format(tag=packVer.nameVer_.version)
                try:
                    Utils.run(cCmd)
                    Utils.run_in_dir(tagCmd, bPath.build_sub_dir)
                except Exception, e:
                    print e
                    print "failed to clone from {url}".format(url=bPath.url)
                    sys.exit(1)
            try:
                Utils.run_in_dir(cmd, bPath.build_sub_dir)
            except Exception, e:
                print e
                print "failed to build in {BUILD}, removing {LOCAL}".format(BUILD=bPath.build_sub_dir, LOCAL = bPath.local_dir)
                Utils.rm_rf(bPath.local_dir)
                sys.exit(1)
    
    def __gitBranch(self, packVer):
        bPath = packVer.bPaths_
        '''
            For header only libraries, will be put directly into local
        '''
        if self.noInternet_:
            newUrl = bPath.url.replace(".git","/archive/" + str(packVer.nameVer_.version) + ".tar.gz").replace("git@github.com:", "https://github.com/")
            base_file = os.path.basename(newUrl)
            fnp = os.path.join(self.dirMaster_.ext_tars,packVer.nameVer_.name, base_file)
            Utils.clear_dir(os.path.dirname(bPath.local_dir))
            Utils.untar(fnp, os.path.dirname(bPath.local_dir))
            ## might not be the best way to do this but works for now
            untaredDir = os.listdir(os.path.dirname(bPath.local_dir))[0]
            os.rename(os.path.join(os.path.dirname(bPath.local_dir), untaredDir), bPath.local_dir)
        else:
            print "cloning from {url}".format(url=bPath.url)
            cCmd = "git clone -b {branch} {url} {d}".format(branch = packVer.nameVer_.version.replace("__", "/"),url=bPath.url, d=bPath.local_dir)
            try:
                Utils.run(cCmd)
            except Exception, e:
                print e
                print "failed to clone branch {branch} from {url}".format(branch = packVer.nameVer_.version.replace("__", "/"), url=bPath.url)
                sys.exit(1)
    
    def __gitTag(self, packVer):
        bPath = packVer.bPaths_
        '''
            For header only libraries, will be put directly into local
        '''
        if self.noInternet_:
            newUrl = bPath.url.replace(".git","/archive/" + str(packVer.nameVer_.version) + ".tar.gz").replace("git@github.com:", "https://github.com/")
            base_file = os.path.basename(newUrl)
            fnp = os.path.join(self.dirMaster_.ext_tars,packVer.nameVer_.name, base_file)
            Utils.clear_dir(os.path.dirname(bPath.local_dir))
            Utils.untar(fnp, os.path.dirname(bPath.local_dir))
            ## might not be the best way to do this but works for now
            untaredDir = os.listdir(os.path.dirname(bPath.local_dir))[0]
            os.rename(os.path.join(os.path.dirname(bPath.local_dir), untaredDir), bPath.local_dir)
        else:
            cmd = "git clone {url} {d}".format(url=bPath.url, d=Utils.shellquote(bPath.local_dir))
            tagCmd = "git checkout {tag}".format(tag=packVer.nameVer_.version)
            try:
                Utils.run(cmd)
                Utils.run_in_dir(tagCmd, bPath.local_dir)
            except:
                print "failed to clone from {url}".format(url=bPath.url)
                sys.exit(1)
    
    def __defaultBuild(self, package, version, fromGitTag = True):
        pack = self.__package(package)
        if not pack.hasVersion(version):
            raise Exception("No set up for version " + str(version) + " for " + str(package))
        packVer = pack.versions_[version]
        bPaths = packVer.bPaths_
        cmd = pack.defaultBuildCmd_.format(external = Utils.shellquote(self.dirMaster_.base_dir), build_sub_dir = Utils.shellquote(bPaths.build_sub_dir), local_dir=Utils.shellquote(bPaths.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        Utils.mkdir(os.path.dirname(bPaths.local_dir))
        if "" != cmd and self.args.verbose:
            print cmd
        if "git" == pack.libType_:
            Utils.mkdir(bPaths.build_dir)
            if fromGitTag:
                self.__buildFromGitTag(packVer, cmd)
            else:
                self.__buildFromGitBranch(packVer, cmd)
        elif "git-headeronly" == pack.libType_:
            if fromGitTag:
                self.__gitTag(packVer)
            else:
                self.__gitBranch(packVer)
        elif "file" == pack.libType_:
            Utils.mkdir(bPaths.build_dir)
            self.__buildFromFile(packVer, cmd)
        elif "file-executable" == pack.libType_:
            Utils.mkdir(bPaths.build_dir)
            self.__buildFromFileExecutable(packVer, cmd)
        else:
            raise Exception("Unrecognized lib type " + str(pack.libType_))
        if Utils.isMac():
            libPath = os.path.join(bPaths.local_dir, "lib")
            if(os.path.exists(libPath)):
                Utils.fixDyLibOnMac(libPath)

            
    def linkInBin(self, package, version, overwrite = False):
        self.packages_.checkForPackVer(LibNameVer(package, version))
        masterBinDir = os.path.join(os.path.dirname(self.extDirLoc), "bin" )
        Utils.mkdir(masterBinDir)
        masterBinDir = os.path.abspath(masterBinDir)
        pack = self.packages_.package(package)
        installDir = pack.getLocalDir(version)
        if os.path.exists(os.path.join(installDir, "bin")):
            binFiles = os.listdir(os.path.join(installDir, "bin"))
            for bFile in binFiles:
                bFileFull = os.path.join(installDir, "bin", bFile)
                if os.path.isfile(bFileFull) and os.access(bFileFull, os.X_OK):
                    if os.path.exists(os.path.join(masterBinDir, bFile)):
                        if overwrite:
                            os.remove(os.path.join(masterBinDir, bFile))
                        else:
                            raise Exception("File: " + os.path.join(masterBinDir, bFile) + " already exists, use --overWrite to overWrite")
                    print "Linking " + CT.boldGreen(bFileFull) + " to " + CT.boldBlue(os.path.join(masterBinDir, bFile))
                    os.symlink(bFileFull, os.path.join(masterBinDir, bFile))



    def armadillo(self, version):
        self.__defaultBuild("armadillo", version)

    
    def downloadFiles(self):
        for set in self.setUpsNeeded:
            topTempDir = os.path.join(self.dirMaster_.base_dir, "temp")
            self.packages_.checkForPackVer(set)
            pack = self.__package(set.name) 
            packVer = pack.versions_[set.version]
            downloadDir = os.path.join(self.dirMaster_.ext_tars, pack.name_)
            Utils.mkdir(downloadDir)
            if pack.bibProject_:
                downloadCmd = "python ./configure.py -CC {CC} -CXX {CXX} -externalLibDir {external} && ./setup.py --compfile compfile.mk --justDownload".format(external = Utils.shellquote(self.dirMaster_.base_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
                Utils.mkdir(topTempDir)
                packVer = pack.versions_[set.version]
                tempDir = os.path.join(topTempDir, pack.name_)
                cloneCmd = "git clone {url} {d}".format(url=packVer.bPaths_.url, d = tempDir)
                tagCmd = "git checkout {tag}".format(tag=packVer.nameVer_.version.replace("__", "/"))
                Utils.run(cloneCmd)
                Utils.run_in_dir(tagCmd, tempDir)
                Utils.run_in_dir(downloadCmd, tempDir)
                if "develop" == set.version or "master" == set.version or "release" in set.version:
                    archiveCmd = "git archive --prefix={name}/ -o {downloadDir}/{version}.tar.gz HEAD".format(name = pack.name_, downloadDir = downloadDir, version = set.version)
                    Utils.run_in_dir(archiveCmd, tempDir)
                shutil.rmtree(tempDir)
            if pack.bibProject_ and ("develop" == set.version or "master" == set.version or "release" in set.version):
                pass
            else:
                url = packVer.getDownloadUrl()
                dest = os.path.join(self.dirMaster_.ext_tars, packVer.nameVer_.name)
                print ("Downloading " + CT.boldGreen(url) + " to " + CT.boldBlue(dest))
                if pack.libType_.startswith("git"):
                    fnp = Utils.get_file(url, dest)
                else:
                    fnp = Utils.get_file_if_size_diff(url, dest)
                
        if os.path.exists(os.path.join(self.dirMaster_.base_dir, "temp")) and os.listdir(os.path.join(self.dirMaster_.base_dir, "temp")) == []:
            shutil.rmtree(os.path.join(self.dirMaster_.base_dir, "temp"))
        print ("Now run \"./setup.py --compfile compfile.mk --outMakefile makefile-common.mk --noInternet\" to build libraries")

    def externalChecks(self):
        ccWhich = Utils.which(self.CC)
        cxxWhich = Utils.which(self.CXX)
        cmakeWhich = Utils.which("cmake")
        gitWhich = Utils.which("git")
        if not ccWhich or not cxxWhich or not cmakeWhich or not gitWhich:
            if not ccWhich:
                print CT.boldRed("Could not find c compiler " + CT.purple + self.CC)
                if self.args.compfile:
                    print "Change CC in " + self.args.compfile
                else:
                    print "Can supply another c compiler by using -CC [option] or by defining bash environmental CC "
                print ""
            if not cxxWhich:
                print CT.boldRed("Could not find c++ compiler " + CT.purple + self.CXX)
                if self.args.compfile:
                    print "Change CXX in " + self.args.compfile
                else:
                    print "Can supply another c++ compiler by using -CXX [option] or by defining bash environmental CXX "
                print ""
            if not cmakeWhich:
                print CT.boldRed("Could not find " + CT.purple + "cmake")
                if Utils.isMac():
                    print "If you have brew, you can install via, brew update && brew install cmake, otherwise you can follow instructions from http://www.cmake.org/install/"
                else:
                    print "On ubuntu to install latest cmake do the following"
                    print "sudo add-apt-repository ppa:george-edison55/cmake-3.x"
                    print "sudo apt-get update"
                    print "sudo apt-get install cmake"
            if not gitWhich:
                print "Can't find git"
            raise Exception("")
        
    def clearCache(self):
        Utils.rm_rf(self.dirMaster_.cache_dir)
        Utils.mkdir(self.dirMaster_.cache_dir)

def ubuntu(self):
        pkgs = """libbz2-dev python2.7-dev cmake libpcre3-dev zlib1g-dev libgcrypt11-dev libicu-dev
python doxygen doxygen-gui auctex xindy graphviz libcurl4-openssl-dev""".split()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compfile', type=str, nargs=1)
    parser.add_argument('--libs', type=str, help="The libraries to install")
    parser.add_argument('--printLibs', action = "store_true", help="Print Available Libs")
    parser.add_argument('--printGitRefs', action = "store_true", help="Print Git branhes and tags for git projects")
    parser.add_argument('--forceUpdate', action = "store_true", help="Remove already installed libs and re-install")
    parser.add_argument('--CC', type=str, nargs=1)
    parser.add_argument('--CXX', type=str, nargs=1)
    parser.add_argument('--numCores', type=str)
    parser.add_argument('--outMakefile', type=str)
    parser.add_argument('--overWrite', action = 'store_true')
    parser.add_argument('--append', action = 'store_true')
    parser.add_argument('--noInternet', action = 'store_true')
    parser.add_argument('--justDownload', action = 'store_true')
    parser.add_argument('--verbose', action = 'store_true')
    parser.add_argument('--symlinkBin', action = 'store_true', help = "Symlink in executables into a directory bin next to external")
    parser.add_argument('--clearCache', action = 'store_true')
    return parser.parse_args()



def main():
    args = parse_args()
    s = Setup(args)
    s.externalChecks()

    if args.printLibs:
        s.printAvailableSetUps()
        return 0
    else:
        if len(s.setUpsNeeded) == 0 and not args.compfile:
            s.printAvailableSetUps()
            return 1
        elif args.printGitRefs:
            s.printGitRefs()
            return 0
        else:
            if args.justDownload:
                s.downloadFiles()
            else:
                s.setup()
                
                if args.outMakefile:
                    packVers = []
                    for set in s.setUpsNeeded:
                        s.packages_.addPackage(packVers,set)
                    s.packages_.writeMakefile(packVers, args.outMakefile, args.overWrite, args.append)
                if args.symlinkBin:
                    for set in s.setUpsNeeded:
                        s.linkInBin(set.name, set.version, args.overWrite)
                return 0

if __name__ == '__main__':
    main()
    
    
    
