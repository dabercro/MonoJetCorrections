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

OutputPerpName      = config.get('OutputBranches','OutputPerpName')
OutputParaName      = config.get('OutputBranches','OutputParaName')
OutputRecoilPtName  = config.get('OutputBranches','OutputRecoilPtName')
OutputParaScaleUpName      = config.get('OutputBranches','OutputParaScaleUpName')
OutputRecoilPtScaleUpName  = config.get('OutputBranches','OutputRecoilPtScaleUpName')
OutputParaScaleDownName      = config.get('OutputBranches','OutputParaScaleDownName')
OutputRecoilPtScaleDownName  = config.get('OutputBranches','OutputRecoilPtScaleDownName')
OutputPerpResolutionUpName      = config.get('OutputBranches','OutputPerpResolutionUpName')
OutputParaResolutionUpName      = config.get('OutputBranches','OutputParaResolutionUpName')
OutputRecoilPtResolutionUpName  = config.get('OutputBranches','OutputRecoilPtResolutionUpName')
OutputPerpResolutionDownName      = config.get('OutputBranches','OutputPerpResolutionDownName')
OutputParaResolutionDownName      = config.get('OutputBranches','OutputParaResolutionDownName')
OutputRecoilPtResolutionDownName  = config.get('OutputBranches','OutputRecoilPtResolutionDownName')
OutputParaFootUpName      = config.get('OutputBranches','OutputParaFootUpName')
OutputRecoilPtFootUpName  = config.get('OutputBranches','OutputRecoilPtFootUpName')
OutputParaFootDownName      = config.get('OutputBranches','OutputParaFootDownName')
OutputRecoilPtFootDownName  = config.get('OutputBranches','OutputRecoilPtFootDownName')
reportFreq  = int(config.get('Misc','ReportFrequency'))

OutputPerpLinearUpName      = config.get('OutputBranches','OutputPerpLinearUpName')
OutputParaLinearUpName      = config.get('OutputBranches','OutputParaLinearUpName')
OutputRecoilPtLinearUpName  = config.get('OutputBranches','OutputRecoilPtLinearUpName')
OutputPerpLinearDownName      = config.get('OutputBranches','OutputPerpLinearDownName')
OutputParaLinearDownName      = config.get('OutputBranches','OutputParaLinearDownName')
OutputRecoilPtLinearDownName  = config.get('OutputBranches','OutputRecoilPtLinearDownName')
OutputPerpQuadUpName      = config.get('OutputBranches','OutputPerpQuadUpName')
OutputParaQuadUpName      = config.get('OutputBranches','OutputParaQuadUpName')
OutputRecoilPtQuadUpName  = config.get('OutputBranches','OutputRecoilPtQuadUpName')
OutputPerpQuadDownName      = config.get('OutputBranches','OutputPerpQuadDownName')
OutputParaQuadDownName      = config.get('OutputBranches','OutputParaQuadDownName')
OutputRecoilPtQuadDownName  = config.get('OutputBranches','OutputRecoilPtQuadDownName')

if config.has_option('InputExpressions','MacrosToLoad'):
    macros = (config.get('InputExpressions','MacrosToLoad')).strip(' ').split(',')
else:
    macros = []
##

for aMacro in macros:
    print "Loading " + aMacro + " ..."
    ROOT.gROOT.LoadMacro('macros/'+aMacro)
##

ROOT.gROOT.LoadMacro('macros/RecoilCorrector.cc+')

if len(macros) == 0:
    ROOT.gROOT.LoadMacro('macros/applicatorOnlyBranches.cc+')
    applicator = ROOT.applicatorOnlyBranches
else:
    ROOT.gROOT.LoadMacro('macros/applicator.cc+')
    applicator = ROOT.applicator
##

sqrt = ROOT.TMath.Sqrt

phoCorrections = ROOT.TFile("data/FootprintFits.root")
ZmmFunc = phoCorrections.Get("mu_Zmm_Data")
ZmmFuncUp = phoCorrections.Get("mu_up_Zmm_Data")
ZmmFuncDown = phoCorrections.Get("mu_down_Zmm_Data")
ZeeFunc = phoCorrections.Get("mu_Zee_Data")
ZeeFuncUp = phoCorrections.Get("mu_up_Zee_Data")
ZeeFuncDown = phoCorrections.Get("mu_down_Zee_Data")
GJetsFunc = phoCorrections.Get("mu_gjets_Data")
GJetsFuncUp = phoCorrections.Get("mu_up_gjets_Data")
GJetsFuncDown = phoCorrections.Get("mu_down_gjets_Data")

def ApplyCorrection(inQueue):
    running = True
    smearingCorrections = ROOT.TFile("data/SmearingFits.root")
    rc = ROOT.RecoilCorrector()
    rc.SetInputName("Zmm")
    rc.LoadAllFits(smearingCorrections)
    while running:
        try:
            inFileName = inQueue.get(True,2)
            print "About to process " + inFileName
            startTime = time()
            applicator(
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
                OutputPerpName,
                OutputParaName,
                OutputRecoilPtName,
                OutputParaScaleUpName,
                OutputRecoilPtScaleUpName,
                OutputParaScaleDownName,
                OutputRecoilPtScaleDownName,
                OutputPerpResolutionUpName,
                OutputParaResolutionUpName,
                OutputRecoilPtResolutionUpName,
                OutputPerpResolutionDownName,
                OutputParaResolutionDownName,
                OutputRecoilPtResolutionDownName,
                OutputParaFootUpName,
                OutputRecoilPtFootUpName,
                OutputParaFootDownName,
                OutputRecoilPtFootDownName,
                OutputPerpLinearUpName,
                OutputParaLinearUpName,
                OutputRecoilPtLinearUpName,
                OutputPerpLinearDownName,
                OutputParaLinearDownName,
                OutputRecoilPtLinearDownName,
                OutputPerpQuadUpName,
                OutputParaQuadUpName,
                OutputRecoilPtQuadUpName,
                OutputPerpQuadDownName,
                OutputParaQuadDownName,
                OutputRecoilPtQuadDownName,
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
