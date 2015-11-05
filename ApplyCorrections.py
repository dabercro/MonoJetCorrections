#! /usr/bin/python

import ROOT
import os, ConfigParser, warnings, sys
from array import array
from multiprocessing import Process, Queue
from optparse import OptionParser
from time import time

warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )

parser = OptionParser()

parser.add_option('-c',help='Names a configuration file.',dest='configName',action='store',metavar='<name>')

(opts,args) = parser.parse_args()

cfgName = "default.cfg"
if opts.__dict__['configName']:
    cfgName = opts.configName

config = ConfigParser.RawConfigParser()
config.read(cfgName)
numMaxProcesses = int(config.get('Misc','NumMaxProcesses'))
inDir           = config.get('FileLocation','InDir')
outDir          = config.get('FileLocation','OutDir')
inTreeName      = config.get('TreeNames','InTreeName')
outTreeName     = config.get('TreeNames','OutTreeName')
PhotonPtExpression      = config.get('InputExpressions','PhotonPtExpression')
DiLeptonPtExpression    = config.get('InputExpressions','DiLeptonPtExpression')
GenBosonPtExpression    = config.get('InputExpressions','GenBosonPtExpression')
GenBosonPdgIdExpression = config.get('InputExpressions','GenBosonPdgIdExpression')
DaughterPdgIdExpression = config.get('InputExpressions','DaughterPdgIdExpression')
EventNumExpression      = config.get('InputExpressions','EventNumExpression')
OutputPhotonPt      = config.get('OutputBranches','OutputPhotonPt')
OutputPhotonSysUp   = config.get('OutputBranches','OutputPhotonSysUp')
OutputPhotonSysDown = config.get('OutputBranches','OutputPhotonSysDown')
OutputZPt      = config.get('OutputBranches','OutputZPt')
OutputZSysUp   = config.get('OutputBranches','OutputZSysUp')
OutputZSysDown = config.get('OutputBranches','OutputZSysDown')
OutputPerpName      = config.get('OutputBranches','OutputPerpName')
OutputParaName      = config.get('OutputBranches','OutputParaName')
OutputRecoilPtName  = config.get('OutputBranches','OutputRecoilPtName')
OutputPerpUpName      = config.get('OutputBranches','OutputPerpUpName')
OutputParaUpName      = config.get('OutputBranches','OutputParaUpName')
OutputRecoilPtUpName  = config.get('OutputBranches','OutputRecoilPtUpName')
OutputPerpDownName      = config.get('OutputBranches','OutputPerpDownName')
OutputParaDownName      = config.get('OutputBranches','OutputParaDownName')
OutputRecoilPtDownName  = config.get('OutputBranches','OutputRecoilPtDownName')
reportFreq  = int(config.get('Misc','ReportFrequency'))

if config.has_option('InputExpressions','MacrosToLoad'):
    macros = (config.get('InputExpressions','MacrosToLoad')).strip(' ').split(',')
else:
    macros = []
##

for aMacro in macros:
    print "Loading " + aMacro + " ..."
    ROOT.gROOT.LoadMacro(aMacro)
##

ROOT.gROOT.LoadMacro('RecoilCorrector.cc+')
ROOT.gROOT.LoadMacro('applicator.cc+')

sqrt = ROOT.TMath.Sqrt

phoCorrections = ROOT.TFile("FootprintFits.root")
ZmmFunc = phoCorrections.Get("mu_Zmm_Data")
ZmmFuncUp = phoCorrections.Get("mu_up_Zmm_Data")
ZmmFuncDown = phoCorrections.Get("mu_down_Zmm_Data")
ZeeFunc = phoCorrections.Get("mu_Zee_Data")
ZeeFuncUp = phoCorrections.Get("mu_up_Zee_Data")
ZeeFuncDown = phoCorrections.Get("mu_down_Zee_Data")
GJetsFunc = phoCorrections.Get("mu_gjets_Data")
GJetsFuncUp = phoCorrections.Get("mu_up_gjets_Data")
GJetsFuncDown = phoCorrections.Get("mu_down_gjets_Data")

branchesCheck = [
    OutputPhotonPt,
    OutputPhotonSysUp,
    OutputPhotonSysDown,
    OutputZPt,
    OutputZSysUp,
    OutputZSysDown,
    OutputPerpName,
    OutputParaName,
    OutputRecoilPtName,
    OutputPerpUpName,
    OutputParaUpName,
    OutputRecoilPtUpName,
    OutputPerpDownName,
    OutputParaDownName,
    OutputRecoilPtDownName,
]

def ApplyCorrection(inQueue):
    running = True
    smearingCorrections = ROOT.TFile("SmearingFits.root")
    rc = ROOT.RecoilCorrector()
    rc.SetInputName("Zmm")
    rc.LoadAllFits(smearingCorrections)
    while running:
        try:
            inFileName = inQueue.get(True,2)
            print "About to process " + inFileName
            startTime = time()
            ROOT.applicator(
                outDir, 
                inDir,
                inFileName,
                inTreeName,
                outTreeName,
                PhotonPtExpression,
                DiLeptonPtExpression,
                GenBosonPtExpression,
                GenBosonPdgIdExpression,
                DaughterPdgIdExpression,
                EventNumExpression,
                OutputPhotonPt,
                OutputPhotonSysUp,
                OutputPhotonSysDown,
                OutputZPt,
                OutputZSysUp,
                OutputZSysDown,
                OutputPerpName,
                OutputParaName,
                OutputRecoilPtName,
                OutputPerpUpName,
                OutputParaUpName,
                OutputRecoilPtUpName,
                OutputPerpDownName,
                OutputParaDownName,
                OutputRecoilPtDownName,
                ZmmFunc,
                ZmmFuncUp,
                ZmmFuncDown,
                ZeeFunc,
                ZeeFuncUp,
                ZeeFuncDown,
                GJetsFunc,
                GJetsFuncUp,
                GJetsFuncDown,
                rc,
                reportFreq
                )
            print "Finished " + inFileName + " ... Elapsed time: " + str(time() - startTime) + " seconds"
            ##
        except:
            print "Worker finished ..."
            smearingCorrections.Close()
            running = False
        ##
    ##
##

theQueue     = Queue()
theProcesses = []

for inFileName in os.listdir(inDir):
    if inFileName.endswith(".root"):
        theQueue.put(inFileName)
##

print ""
print "About to start " + str(numMaxProcesses) + " workers!"
print ""

for worker in range(numMaxProcesses):
    aProcess = Process(target=ApplyCorrection, args=(theQueue,))
    aProcess.start()
    theProcesses.append(aProcess)
##

for aProccess in theProcesses:
    aProccess.join()
##

phoCorrections.Close()

print "All done!"
