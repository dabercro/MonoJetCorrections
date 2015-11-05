#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"
#include "TMath.h"

#include "RecoilCorrector.h"

void applicatorAllBranches(TString outDir, 
                           TString inDir,
                           TString inFileName,
                           TString inTreeName,
                           TString outTreeName,
                           TString PhotonPtExpression,
                           TString DiLeptonPtExpression,
                           TString GenBosonPtExpression,
                           TString GenBosonPdgIdExpression,
                           TString DaughterPdgIdExpression,
                           TString EventNumExpression,
                           TString OutputPhotonPt,
                           TString OutputPhotonSysUp,
                           TString OutputPhotonSysDown,
                           TString OutputZPt,
                           TString OutputZSysUp,
                           TString OutputZSysDown,
                           TString OutputPerpName,
                           TString OutputParaName,
                           TString OutputRecoilPtName,
                           TString OutputPerpUpName,
                           TString OutputParaUpName,
                           TString OutputRecoilPtUpName,
                           TString OutputPerpDownName,
                           TString OutputParaDownName,
                           TString OutputRecoilPtDownName,
                           TF1* ZmmFunc,
                           TF1* ZmmFuncUp,
                           TF1* ZmmFuncDown,
                           TF1* ZeeFunc,
                           TF1* ZeeFuncUp,
                           TF1* ZeeFuncDown,
                           TF1* GJetsFunc,
                           TF1* GJetsFuncUp,
                           TF1* GJetsFuncDown,
                           RecoilCorrector* rc,
                           int reportFreq
                           )

{
  TFile *inFile  = 0;
  TFile *outFile = 0;
  if (outDir == inDir) {
    inFile  = new TFile(inDir + TString("/") + inFileName,"UPDATE");
    outFile = inFile;
  }
  else {
    inFile  = new TFile(inDir + TString("/") + inFileName);
    outFile = new TFile(outDir + TString("/") + inFileName,"RECREATE");
  }
  TTree *inTree  = (TTree*) inFile->Get(inTreeName);
  if (inTree == NULL) {
    std::cout <<  "" << std::endl;
    std::cout <<  "###############################################################" << std::endl;
    std::cout <<  "##" << std::endl;
    std::cout <<  "##   InTree not successfully fetched in " << inFileName << "!" << std::endl;
    std::cout <<  "##   Check config file for tree name." << std::endl;
    std::cout <<  "##" << std::endl;
    std::cout <<  "##   Aborting ..." << std::endl;
    std::cout <<  "##" << std::endl;
    std::cout <<  "###############################################################" << std::endl;
    std::cout <<  "" << std::endl;
    outFile->Close();
    if (inFile->IsOpen())
      inFile->Close();
    return;
  }
  TTree *outTree = 0;
  if (outTreeName == inTreeName && outFile == inFile)
    outTree = inTree;
  else {
    outFile->cd();
    outTree = new TTree(outTreeName,outTreeName);
  }

  float photonPt = 0;
  float ZPt      = 0;

  float genBosPt    = 0;
  int genBosPdgId   = 0;
  int daughterPdgId = 0;
  int eventNum      = 0;

  TBranch *photonPtF = inTree->GetBranch(PhotonPtExpression);
  photonPtF->SetAddress(&photonPt);
  TBranch *ZPtF      = inTree->GetBranch(DiLeptonPtExpression);
  ZPtF->SetAddress(&ZPt);

  TBranch *genBosPtF      = inTree->GetBranch(GenBosonPtExpression);
  genBosPtF->SetAddress(&genBosPt);
  TBranch *genBosPdgIdF   = inTree->GetBranch(GenBosonPdgIdExpression);
  genBosPdgIdF->SetAddress(&genBosPdgId);
  TBranch *daughterPdgIdF = inTree->GetBranch(DaughterPdgIdExpression);
  daughterPdgIdF->SetAddress(&daughterPdgId);
  TBranch *eventNumF      = inTree->GetBranch(EventNumExpression);
  eventNumF->SetAddress(&eventNum);

  float footPt = 0;
  float footUp = 0;
  float footDown = 0;
  float ZfootPt = 0;
  float ZfootUp = 0;
  float ZfootDown = 0;
  float uPerp = 0;
  float uPara = 0;
  float uMag = 0;
  float uPerpUp = 0;
  float uParaUp = 0;
  float uMagUp = 0;
  float uPerpDown = 0;
  float uParaDown = 0;
  float uMagDown = 0;

  TBranch *footPtBr    = outTree->Branch(OutputPhotonPt,&footPt,OutputPhotonPt+"/F");
  TBranch *footUpBr    = outTree->Branch(OutputPhotonSysUp,&footUp,OutputPhotonSysUp+"/F");
  TBranch *footDownBr  = outTree->Branch(OutputPhotonSysDown,&footDown,OutputPhotonSysDown+"/F");
  TBranch *ZfootPtBr   = outTree->Branch(OutputZPt,&ZfootPt,OutputZPt+"/F");
  TBranch *ZfootUpBr   = outTree->Branch(OutputZSysUp,&ZfootUp,OutputZSysUp+"/F");
  TBranch *ZfootDownBr = outTree->Branch(OutputZSysDown,&ZfootDown,OutputZSysDown+"/F");

  TBranch *uPerpBr = outTree->Branch(OutputPerpName,&uPerp,OutputPerpName+"/F");
  TBranch *uParaBr = outTree->Branch(OutputParaName,&uPara,OutputParaName+"/F");
  TBranch *uMagBr  = outTree->Branch(OutputRecoilPtName,&uMag,OutputRecoilPtName+"/F");
  TBranch *uPerpUpBr = outTree->Branch(OutputPerpUpName,&uPerpUp,OutputPerpUpName+"/F");
  TBranch *uParaUpBr = outTree->Branch(OutputParaUpName,&uParaUp,OutputParaUpName+"/F");
  TBranch *uMagUpBr  = outTree->Branch(OutputRecoilPtUpName,&uMagUp,OutputRecoilPtUpName+"/F");
  TBranch *uPerpDownBr = outTree->Branch(OutputPerpDownName,&uPerpDown,OutputPerpDownName+"/F");
  TBranch *uParaDownBr = outTree->Branch(OutputParaDownName,&uParaDown,OutputParaDownName+"/F");
  TBranch *uMagDownBr  = outTree->Branch(OutputRecoilPtDownName,&uMagDown,OutputRecoilPtDownName+"/F");

  int lastPdgId = 0;

  int numEntries = inTree->GetEntriesFast();
  for (int entry = 0; entry != numEntries; ++entry) {
    if (entry % reportFreq == 0)
      std::cout <<  "Processing " << inFileName << " ... " << float(entry)/numEntries * 100 << "%" << std::endl;

    photonPtF->GetEntry(entry);
    ZPtF->GetEntry(entry);
    genBosPtF->GetEntry(entry);
    genBosPdgIdF->GetEntry(entry);
    daughterPdgIdF->GetEntry(entry);
    eventNumF->GetEntry(entry);

    if (photonPt > 1000 || photonPt < 0) {
        footPt   = photonPt;
        footUp   = photonPt;
        footDown = photonPt;
      }
    else {
      double mmFuncVal  = ZmmFunc->Eval(photonPt);
      double mmUpErr2   = ZmmFuncUp->Eval(photonPt) - ZmmFunc->Eval(photonPt);
      double mmDownErr2 = ZmmFunc->Eval(photonPt) - ZmmFuncDown->Eval(photonPt);

      double gjetsFuncVal  = GJetsFunc->Eval(photonPt);
      double gjetsUpErr2   = GJetsFuncUp->Eval(photonPt) - GJetsFunc->Eval(photonPt);
      double gjetsDownErr2 = GJetsFunc->Eval(photonPt) - GJetsFuncDown->Eval(photonPt);

      footPt   = photonPt + (mmFuncVal - gjetsFuncVal)/(1 - ZmmFunc->GetParameter(1));
      footUp   = photonPt + (mmFuncVal - gjetsFuncVal + sqrt(mmUpErr2 + gjetsDownErr2))/(1 - ZmmFuncUp->GetParameter(1) - 2*ZmmFuncUp->GetParameter(2)*photonPt);
      footDown = photonPt + (mmFuncVal - gjetsFuncVal - sqrt(mmDownErr2 + gjetsUpErr2))/(1 - ZmmFuncDown->GetParameter(1) - 2*ZmmFuncDown->GetParameter(2)*photonPt);
    }
    if (ZPt > 1000 || ZPt < 0 || abs(daughterPdgId) != 11) {
      ZfootPt   = ZPt;
      ZfootUp   = ZPt;
      ZfootDown = ZPt;
    }
    else {
      double mmFuncVal  = ZmmFunc->Eval(photonPt);
      double mmUpErr2   = ZmmFuncUp->Eval(photonPt) - ZmmFunc->Eval(photonPt);
      double mmDownErr2 = ZmmFunc->Eval(photonPt) - ZmmFuncDown->Eval(photonPt);

      double eeFuncVal  = ZeeFunc->Eval(photonPt);
      double eeUpErr2   = ZeeFuncUp->Eval(photonPt) - ZeeFunc->Eval(photonPt);
      double eeDownErr2 = ZeeFunc->Eval(photonPt) - ZeeFuncDown->Eval(photonPt);

      footPt   = ZPt + (mmFuncVal - eeFuncVal)/(1 - ZmmFunc->GetParameter(1));
      footUp   = ZPt + (mmFuncVal - eeFuncVal + sqrt(mmUpErr2 + eeDownErr2))/(1 - ZmmFuncUp->GetParameter(1) - 2*ZmmFuncUp->GetParameter(2)*ZPt);
      footDown = ZPt + (mmFuncVal - eeFuncVal - sqrt(mmDownErr2 + eeUpErr2))/(1 - ZmmFuncDown->GetParameter(1) - 2*ZmmFuncDown->GetParameter(2)*ZPt);
    }

    rc->SetSeed(eventNum);

    if (genBosPdgId == 22 || genBosPdgId == 23 || abs(genBosPdgId) == 24) {
      if (genBosPdgId != lastPdgId) {
        lastPdgId = genBosPdgId;
        if (genBosPdgId == 22)
          rc->SetOutput(RecoilCorrector::kGJets);
        else if (genBosPdgId == 23) {
          if (abs(daughterPdgId) == 11)
            rc->SetOutput(RecoilCorrector::kZee);
          else if (abs(daughterPdgId) == 13)
            rc->SetOutput(RecoilCorrector::kZmm);
          else
            rc->SetOutput(RecoilCorrector::kZnn);
        }
        else {
          if (abs(daughterPdgId) == 11)
            rc->SetOutput(RecoilCorrector::kWen);
          else 
            rc->SetOutput(RecoilCorrector::kWmn);
        }
      }
      rc->ComputeU(genBosPt,uPerp,uPara);
      rc->ComputeU(genBosPt,uPerpUp,uParaUp,1.0);
      rc->ComputeU(genBosPt,uPerpDown,uParaDown,-1.0);
      uMag = TMath::Sqrt(uPerp*uPerp + uPara*uPara);
      uMagUp = TMath::Sqrt(uPerpUp*uPerpUp + uParaUp*uParaUp);
      uMagDown = TMath::Sqrt(uPerpDown*uPerpDown + uParaDown*uParaDown);
    }
    else {
      uPerp = 300;
      uPara = 300;
      uMag  = -5;
      uPerpUp = 300;
      uParaUp = 300;
      uMagUp  = -5;
      uPerpDown = 300;
      uParaDown = 300;
      uMagDown  = -5;
    }

    if (inTree == outTree) {
      footPtBr->Fill();
      footUpBr->Fill();
      footDownBr->Fill();
      ZfootPtBr->Fill();
      ZfootUpBr->Fill();
      ZfootDownBr->Fill();
      uPerpBr->Fill();
      uParaBr->Fill();
      uMagBr->Fill();
      uPerpUpBr->Fill();
      uParaUpBr->Fill();
      uMagUpBr->Fill();
      uPerpDownBr->Fill();
      uParaDownBr->Fill();
      uMagDownBr->Fill();
    }
    else
      outTree->Fill();
  }

  outFile->cd();
  outTree->Write();
  outFile->Close();
  if (inFile->IsOpen())
    inFile->Close();
}
