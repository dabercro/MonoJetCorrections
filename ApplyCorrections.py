#! /usr/bin/python

import ROOT
from multiprocessing import Process, Queue
import os
from optparse import OptionParser
import ConfigParser

parser = OptionParser()

parser.add_option('-c',help='Names a configuration file.',dest='configName',action='store',metavar='<name>')

(opts,args) = parser.parse_args()

cfgName = "default.cfg"
if opts.__dict__['configName']:
    cfgName = opts.configName

config = ConfigParser.RawConfigParser()
config.read(cfgName)
numMaxProcesses = config.get('General','NumMaxProcesses')
inDir           = config.get('FileLocation','InDir')
outDir          = config.get('FileLocation','OutDir')
inTreeName      = config.get('TreeNames','InTreeName')
outTreeName     = config.get('TreeNames','OutTreeName')
PhotonPtExpression      = config.get('InputExpressions','PhotonPtExpression')
GenBosonPtExpression    = config.get('InputExpressions','GenBosonPtExpression')
GenBosonPdgIdExpression = config.get('InputExpressions','GenBosonPdgIdExpression')
DaughterPdgIdExpression = config.get('InputExpressions','DaughterPdgIdExpression')
EventNumExpression      = config.get('InputExpressions','EventNumExpression')
OutputPhotonPt       = config.get('OutputBranches','OutputPhotonPt')
OutputPhotonSysUp    = config.get('OutputBranches','OutputPhotonSysUp')
OutputPhotonSysDown  = config.get('OutputBranches','OutputPhotonSysDown')
OutputPerpName       = config.get('OutputBranches','OutputPerpName')
OutputParaName       = config.get('OutputBranches','OutputParaName')
OutputRecoilPtName   = config.get('OutputBranches','OutputRecoilPtName')

ROOT.gROOT.LoadMacro('RecoilCorrector.cc+')

phoCorrections = ROOT.TFile("FootprintFits.root")
ZmmFunc = phoCorrections.Get("mu_Zmm_Data")
GJetsFunc = phoCorrections.Get("mu_gjets_Data")
ZmmFuncUp = phoCorrections.Get("mu_up_Zmm_Data")
GJetsFuncUp = phoCorrections.Get("mu_up_gjets_Data")
ZmmFuncDown = phoCorrections.Get("mu_down_Zmm_Data")
GJetsFuncDown = phoCorrections.Get("mu_down_gjets_Data")

smearingCorrections = ROOT.TFile("SmearingFits.root")

sqrt = ROOT.TMath.Sqrt()

def ApplyCorrection(inQueue):
    running = True
    rc = ROOT.RecoilCorrector()
    rc.SetInName("Zmm")
    while running:
        try:
            inFileIndex = inQueue.get(True,2)
            print "About to process " + inFileName
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
            outTree = 0
            if outTreeName == inTreeName and outFile == inFile:
                outTree = inTree
            else:
                outFile.cd()
                outTree = ROOT.Tree(outTreeName,outTreeName)
            ##
            photonPtF = ROOT.TTreeFormula("PhotonPt",PhotonPtExpression,inTree)

            footPt   = array('f',[0.0])
            footUp   = array('f',[0.0])
            footDown = array('f',[0.0])

            footPtBr   = outTree.Branch(OutputPhotonPt,footPt,OutputPhotonPt+"/F")
            footUpBr   = outTree.Branch(OutputPhotonSysUp,footUp,OutputPhotonSysUp+"/F")
            footDownBr = outTree.Branch(OutputPhotonSysDown,footDown,OutputPhotonSysDown+"/F")

            genBosPtF      = ROOT.TTreeFormula("GenBosPt",GenBosonPtExpression,inTree)
            genBosPdgIdF   = ROOT.TTreeFormula("GenBosPdgId",GenBosonPdgIdExpression,inTree)
            daughterPdgIdF = ROOT.TTreeFormula("DaughterPdgId",DaughterPdgIdExpression,inTree)
            eventNumF      = ROOT.TTreeFormula("EventNum",EventsNumExpression,inTree)

            uPerp = array('f',[0.0])
            uPara = array('f',[0.0])
            uMag  = array('f',[0.0])

            uPerpBr = outTree.Branch(OutputPerpName,uPerp,OutputPerpName+"/F")
            uParaBr = outTree.Branch(OutputParaName,uPara,OutputParaName+"/F")
            uMagBr  = outTree.Branch(OutputRecoilPtName,uMag,OutputRecoilPtName+"/F")

            lastPdgId = 0

            for entry in range(inTree.GetEntriesFast()):
                inTree.GetEntry(entry)

                photonPt = photonPtF.EvalInstance()

                if photonPt > 1000 or photonPt < 0:
                    footPt[0]   = photonPt
                    footUp[0]   = photonPt
                    footDown[0] = photonPt
                ##

                footPt[0]   = photonPt + (ZmmFunc.Eval(photonPt) - GJetsFunc.Eval(photonPt))/(1 - ZmmFunc.GetParameter(1))
                footUp[0]   = photonPt + (ZmmFuncUp.Eval(photonPt) - GJetsFuncDown.Eval(photonPt))/(1 - ZmmFuncUp.GetParameter(1))
                footDown[0] = photonPt + (ZmmFuncDown.Eval(photonPt) - GJetsFuncUp.Eval(photonPt))/(1 - ZmmFuncDown.GetParameter(1))

                genBosPt      = genBosPtF.EvalInstance()
                genBosPdgId   = genBosPdgIdF.EvalInstance()
                daughterPdgId = daughterPdgIdF.EvalInstance()
                eventNum      = eventNumF.EvalInstance()

                rc.SetSeed(eventNum)

                if genBosPdgId in [22,23,24,-24]:
                    if genBosPdgId != lastPdgId:
                        lastPdgId = genBosPdgId
                        if genBosPdg == 22:
                            rc.SetOutputName("GJets")
                        elif genBosPdg == 23:
                            if abs(daughterPdgId) == 11:
                                rc.SetOutputName("Zee")
                            else:
                                rc.SetOutputName("Zee")
                            ##
                        else:
                            if abs(daughterPdgId) == 11:
                                rc.SetOutputName("Wen")
                            else:
                                rc.SetOutputName("Wmn")
                            ##
                        ##
                        rc.LoadAllFits(smearingCorrections)
                    ##
                    rc.ComputeU(genBosPt,uPerp,uPara)
                    uMag[0] = sqrt(uPerp[0]*uPerp[0] + uPara[0] + uPara[0])
                else:
                    uPerp[0] = 300
                    uPara[0] = 300
                    uMag[0]  = -5   ## I'd recommend cutting on this one
                ##
                uPerpBr.Fill()
                uParaBr.Fill()
                uMagBr.Fill()

                footPtBr.Fill()
                footUpBr.Fill()
                footDownBr.Fill()
            ##
            outFile.WriteTObject(outTree,outTreeName)
            outFile.Close()
            if inFile.IsOpen():
                inFile.Close()
            ##
            print "Finished " + inFileName
            ##
        except:
            print "Worker finished..."
            running = False
        ##
    ##
##

theQueue     = Queue()
theProcesses = []

for inFileName in os.listdir(inDir):
    if inFileName.endswith(".root") and not "OLD" in inFileName:
        theQueue.put(inFileName)
##

for worker in range(numMaxProcesses):
    aProcess = Process(target=ApplyCorrection, args=(theQueue,))
    aProcess.start()
    theProcesses.append(aProcess)
##

for aProccess in theProcesses:
    aProccess.join()
##

print "All done!"
