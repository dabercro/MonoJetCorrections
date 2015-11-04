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
            inFile  = 0
            outFile = 0
            if outDir == inDir:
                inFile  = ROOT.TFile(inDir + "/" + inFileName,"UPDATE")
                outFile = inFile
            ##
            else:
                inFile  = ROOT.TFile(inDir + "/" + inFileName)
                outFile = ROOT.TFile(outDir + "/" + inFileName,"RECREATE")
            ##
            inTree  = inFile.Get(inTreeName)
            if inTree == None:
                print ""
                print "###############################################################"
                print "##"
                print "##   InTree not successfully fetched in " + inFileName + "!"
                print "##   Check config file for tree name."
                print "##"
                print "##   Aborting ..."
                print "##"
                print "###############################################################"
                print ""
                outFile.Close()
                if inFile.IsOpen():
                    inFile.Close()
                ##
                exit(1)
            ##
            outTree = 0
            if outTreeName == inTreeName and outFile == inFile:
                outTree = inTree
            else:
                outFile.cd()
                outTree = ROOT.TTree(outTreeName,outTreeName)
            ##

            for check in branchesCheck:
                if outTree.GetBranch(check) != None:
                    print "Branch " + check + " already exists in " + outTreeName + "!" 
                    print "Exiting ..."
                    exit(1)
                ##
            ##

            photonPtF = ROOT.TTreeFormula("PhotonPt",PhotonPtExpression,inTree)
            ZPtF      = ROOT.TTreeFormula("ZPt",DiLeptonPtExpression,inTree)

            footPt    = array('f',[0.0])
            footUp    = array('f',[0.0])
            footDown  = array('f',[0.0])
            ZfootPt   = array('f',[0.0])
            ZfootUp   = array('f',[0.0])
            ZfootDown = array('f',[0.0])

            genBosPtF      = ROOT.TTreeFormula("GenBosPt",GenBosonPtExpression,inTree)
            genBosPdgIdF   = ROOT.TTreeFormula("GenBosPdgId",GenBosonPdgIdExpression,inTree)
            daughterPdgIdF = ROOT.TTreeFormula("DaughterPdgId",DaughterPdgIdExpression,inTree)
            eventNumF      = ROOT.TTreeFormula("EventNum",EventNumExpression,inTree)

            uPerp = array('f',[0.0])
            uPara = array('f',[0.0])
            uMag  = array('f',[0.0])
            uPerpUp = array('f',[0.0])
            uParaUp = array('f',[0.0])
            uMagUp  = array('f',[0.0])
            uPerpDown = array('f',[0.0])
            uParaDown = array('f',[0.0])
            uMagDown  = array('f',[0.0])

            footPtBr    = outTree.Branch(OutputPhotonPt,footPt,OutputPhotonPt+"/F")
            footUpBr    = outTree.Branch(OutputPhotonSysUp,footUp,OutputPhotonSysUp+"/F")
            footDownBr  = outTree.Branch(OutputPhotonSysDown,footDown,OutputPhotonSysDown+"/F")
            ZfootPtBr   = outTree.Branch(OutputZPt,ZfootPt,OutputZPt+"/F")
            ZfootUpBr   = outTree.Branch(OutputZSysUp,ZfootUp,OutputZSysUp+"/F")
            ZfootDownBr = outTree.Branch(OutputZSysDown,ZfootDown,OutputZSysDown+"/F")

            uPerpBr = outTree.Branch(OutputPerpName,uPerp,OutputPerpName+"/F")
            uParaBr = outTree.Branch(OutputParaName,uPara,OutputParaName+"/F")
            uMagBr  = outTree.Branch(OutputRecoilPtName,uMag,OutputRecoilPtName+"/F")
            uPerpUpBr = outTree.Branch(OutputPerpUpName,uPerpUp,OutputPerpUpName+"/F")
            uParaUpBr = outTree.Branch(OutputParaUpName,uParaUp,OutputParaUpName+"/F")
            uMagUpBr  = outTree.Branch(OutputRecoilPtUpName,uMagUp,OutputRecoilPtUpName+"/F")
            uPerpDownBr = outTree.Branch(OutputPerpDownName,uPerpDown,OutputPerpDownName+"/F")
            uParaDownBr = outTree.Branch(OutputParaDownName,uParaDown,OutputParaDownName+"/F")
            uMagDownBr  = outTree.Branch(OutputRecoilPtDownName,uMagDown,OutputRecoilPtDownName+"/F")

            lastPdgId = 0

            numEntries = inTree.GetEntriesFast()
            for entry in range(numEntries):
                if entry % reportFreq == 0:
                    print "Processing " + inFileName + " ... " + str(float(entry)/numEntries * 100) + "%"
                ##
                inTree.GetEntry(entry)

                photonPt = float(photonPtF.EvalInstance())
                ZPt      = float(ZPtF.EvalInstance())

                genBosPt      = float(genBosPtF.EvalInstance())
                genBosPdgId   = int(genBosPdgIdF.EvalInstance())
                daughterPdgId = int(daughterPdgIdF.EvalInstance())
                eventNum      = int(eventNumF.EvalInstance())

                if photonPt > 1000 or photonPt < 0:
                    footPt[0]   = photonPt
                    footUp[0]   = photonPt
                    footDown[0] = photonPt
                else:
                    footPt[0]   = photonPt + (ZmmFunc.Eval(photonPt) - GJetsFunc.Eval(photonPt))/(1 - ZmmFunc.GetParameter(1))
                    footUp[0]   = photonPt + (ZmmFuncUp.Eval(photonPt) - GJetsFuncDown.Eval(photonPt))/(1 - ZmmFuncUp.GetParameter(1))
                    footDown[0] = photonPt + (ZmmFuncDown.Eval(photonPt) - GJetsFuncUp.Eval(photonPt))/(1 - ZmmFuncDown.GetParameter(1))
                ##
                if ZPt > 1000 or ZPt < 0 or abs(daughterPdgId) != 11:
                    ZfootPt[0]   = ZPt
                    ZfootUp[0]   = ZPt
                    ZfootDown[0] = ZPt
                else:
                    ZfootPt[0]   = ZPt + (ZmmFunc.Eval(ZPt) - ZeeFunc.Eval(ZPt))/(1 - ZmmFunc.GetParameter(1))
                    ZfootUp[0]   = ZPt + (ZmmFuncUp.Eval(ZPt) - ZeeFuncDown.Eval(ZPt))/(1 - ZmmFuncUp.GetParameter(1))
                    ZfootDown[0] = ZPt + (ZmmFuncDown.Eval(ZPt) - ZeeFuncUp.Eval(ZPt))/(1 - ZmmFuncDown.GetParameter(1))
                ##

                rc.SetSeed(eventNum)

                if genBosPdgId in [22,23,24,-24]:
                    if genBosPdgId != lastPdgId:
                        lastPdgId = genBosPdgId
                        if genBosPdgId == 22:
                            rc.SetOutput(rc.kGJets)
                        elif genBosPdgId == 23:
                            if abs(daughterPdgId) == 11:
                                rc.SetOutput(rc.kZee)
                            elif abs(daughterPdgId == 13):
                                rc.SetOutput(rc.kZmm)
                            else:
                                rc.SetOutput(rc.kZnn)
                            ##
                        else:
                            if abs(daughterPdgId) == 11:
                                rc.SetOutput(rc.kWen)
                            else:
                                rc.SetOutput(rc.kWmn)
                            ##
                        ##
                    ##
                    rc.ComputeU(genBosPt,uPerp,uPara)
                    rc.ComputeU(genBosPt,uPerpUp,uParaUp,1.0)
                    rc.ComputeU(genBosPt,uPerpDown,uParaDown,-1.0)
                    uMag[0] = sqrt(uPerp[0]*uPerp[0] + uPara[0]*uPara[0])
                    uMagUp[0] = sqrt(uPerpUp[0]*uPerpUp[0] + uParaUp[0]*uParaUp[0])
                    uMagDown[0] = sqrt(uPerpDown[0]*uPerpDown[0] + uParaDown[0]*uParaDown[0])
                else:
                    uPerp[0] = 300
                    uPara[0] = 300
                    uMag[0]  = -5
                    uPerpUp[0] = 300
                    uParaUp[0] = 300
                    uMagUp[0]  = -5
                    uPerpDown[0] = 300
                    uParaDown[0] = 300
                    uMagDown[0]  = -5
                ##

                if inTree == outTree:
                    footPtBr.Fill()
                    footUpBr.Fill()
                    footDownBr.Fill()
                    ZfootPtBr.Fill()
                    ZfootUpBr.Fill()
                    ZfootDownBr.Fill()
                    uPerpBr.Fill()
                    uParaBr.Fill()
                    uMagBr.Fill()
                    uPerpUpBr.Fill()
                    uParaUpBr.Fill()
                    uMagUpBr.Fill()
                    uPerpDownBr.Fill()
                    uParaDownBr.Fill()
                    uMagDownBr.Fill()
                else:
                    outTree.Fill()
                ##
            ##
            outFile.cd()
            outTree.Write()
            outFile.Close()
            if inFile.IsOpen():
                inFile.Close()
            ##
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
