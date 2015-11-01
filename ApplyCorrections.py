#! /usr/bin/python

import ROOT
from multiprocessing import Process, Queue
import os
from Corrector_cff import *

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

            genBosPtF    = ROOT.TTreeFormula("GenBosPt",GenBosonPtExpression,inTree)
            genBosPdgIdF = ROOT.TTreeFormula("GenBosPdgId",GenBosonPdgIdExpression,inTree)
            eventNumF    = ROOT.TTreeFormula("EventNum",EventsNumExpression,inTree)

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

                footPt[0]   = photonPt + (ZmmFunc.Eval(photonPt) - GJetsFunc.Eval(photonPt))/(1 - ZmmFunc.GetParameter(1))
                footUp[0]   = photonPt + (ZmmFuncUp.Eval(photonPt) - GJetsFuncDown.Eval(photonPt))/(1 - ZmmFuncUp.GetParameter(1))
                footDown[0] = photonPt + (ZmmFuncDown.Eval(photonPt) - GJetsFuncUp.Eval(photonPt))/(1 - ZmmFuncDown.GetParameter(1))

                genBosPt    = genBosPtF.EvalInstance()
                genBosPdgId = genBosPdgIdF.EvalInstance()
                eventNum    = eventNumF.EvalInstance()

                rc.SetSeed(eventNum)

                if genBosPdgId in [22,23,24,-24]:
                    if genBosPdgId != lastPdgId:
                        lastPdgId = genBosPdgId
                        if genBosPdg == 22:
                            rc.SetOutputName("GJets")
                        elif genBosPdg == 23:
                            rc.SetOutputName("Zmm")      ### I need to come up with a better thing for this!!! Uh oh!!!
                        else:
                            rc.SetOutputName("Wmn")
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
