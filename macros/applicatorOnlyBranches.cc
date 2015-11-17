#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"
#include "TMath.h"

#include "RecoilCorrector.h"

void applicatorOnlyBranches(TString outDir, 
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
                            TString OutputPerpName,
                            TString OutputParaName,
                            TString OutputRecoilPtName,
                            TString OutputParaScaleUpName,
                            TString OutputRecoilPtScaleUpName,
                            TString OutputParaScaleDownName,
                            TString OutputRecoilPtScaleDownName,
                            TString OutputPerpResolutionUpName,
                            TString OutputParaResolutionUpName,
                            TString OutputRecoilPtResolutionUpName,
                            TString OutputPerpResolutionDownName,
                            TString OutputParaResolutionDownName,
                            TString OutputRecoilPtResolutionDownName,
                            TString OutputParaFootUpName,
                            TString OutputRecoilPtFootUpName,
                            TString OutputParaFootDownName,
                            TString OutputRecoilPtFootDownName,
                            TString OutputPerpLinearUpName,
                            TString OutputParaLinearUpName,
                            TString OutputRecoilPtLinearUpName,
                            TString OutputPerpLinearDownName,
                            TString OutputParaLinearDownName,
                            TString OutputRecoilPtLinearDownName,
                            TString OutputPerpQuadUpName,
                            TString OutputParaQuadUpName,
                            TString OutputRecoilPtQuadUpName,
                            TString OutputPerpQuadDownName,
                            TString OutputParaQuadDownName,
                            TString OutputRecoilPtQuadDownName,
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

  float uPerp = 0;
  float uPara = 0;
  float uMag = 0;
  float uParaScaleUp = 0;
  float uMagScaleUp = 0;
  float uParaScaleDown = 0;
  float uMagScaleDown = 0;
  float uPerpResolutionUp = 0;
  float uParaResolutionUp = 0;
  float uMagResolutionUp = 0;
  float uPerpResolutionDown = 0;
  float uParaResolutionDown = 0;
  float uMagResolutionDown = 0;
  float uParaFootUp = 0;
  float uMagFootUp = 0;
  float uParaFootDown = 0;
  float uMagFootDown = 0;
  float uPerpLinearUp = 0;
  float uParaLinearUp = 0;
  float uMagLinearUp = 0;
  float uPerpLinearDown = 0;
  float uParaLinearDown = 0;
  float uMagLinearDown = 0;
  float uPerpQuadUp = 0;
  float uParaQuadUp = 0;
  float uMagQuadUp = 0;
  float uPerpQuadDown = 0;
  float uParaQuadDown = 0;
  float uMagQuadDown = 0;

  TBranch *uPerpBr = outTree->Branch(OutputPerpName,&uPerp,OutputPerpName+"/F");
  TBranch *uParaBr = outTree->Branch(OutputParaName,&uPara,OutputParaName+"/F");
  TBranch *uMagBr  = outTree->Branch(OutputRecoilPtName,&uMag,OutputRecoilPtName+"/F");
  TBranch *uParaScaleUpBr = outTree->Branch(OutputParaScaleUpName,&uParaScaleUp,OutputParaScaleUpName+"/F");
  TBranch *uMagScaleUpBr  = outTree->Branch(OutputRecoilPtScaleUpName,&uMagScaleUp,OutputRecoilPtScaleUpName+"/F");
  TBranch *uParaScaleDownBr = outTree->Branch(OutputParaScaleDownName,&uParaScaleDown,OutputParaScaleDownName+"/F");
  TBranch *uMagScaleDownBr  = outTree->Branch(OutputRecoilPtScaleDownName,&uMagScaleDown,OutputRecoilPtScaleDownName+"/F");
  TBranch *uPerpResolutionUpBr = outTree->Branch(OutputPerpResolutionUpName,&uPerpResolutionUp,OutputPerpResolutionUpName+"/F");
  TBranch *uParaResolutionUpBr = outTree->Branch(OutputParaResolutionUpName,&uParaResolutionUp,OutputParaResolutionUpName+"/F");
  TBranch *uMagResolutionUpBr  = outTree->Branch(OutputRecoilPtResolutionUpName,&uMagResolutionUp,OutputRecoilPtResolutionUpName+"/F");
  TBranch *uPerpResolutionDownBr = outTree->Branch(OutputPerpResolutionDownName,&uPerpResolutionDown,OutputPerpResolutionDownName+"/F");
  TBranch *uParaResolutionDownBr = outTree->Branch(OutputParaResolutionDownName,&uParaResolutionDown,OutputParaResolutionDownName+"/F");
  TBranch *uMagResolutionDownBr  = outTree->Branch(OutputRecoilPtResolutionDownName,&uMagResolutionDown,OutputRecoilPtResolutionDownName+"/F");
  TBranch *uParaFootUpBr = outTree->Branch(OutputParaFootUpName,&uParaFootUp,OutputParaFootUpName+"/F");
  TBranch *uMagFootUpBr  = outTree->Branch(OutputRecoilPtFootUpName,&uMagFootUp,OutputRecoilPtFootUpName+"/F");
  TBranch *uParaFootDownBr = outTree->Branch(OutputParaFootDownName,&uParaFootDown,OutputParaFootDownName+"/F");
  TBranch *uMagFootDownBr  = outTree->Branch(OutputRecoilPtFootDownName,&uMagFootDown,OutputRecoilPtFootDownName+"/F");
  TBranch *uPerpLinearUpBr = outTree->Branch(OutputPerpLinearUpName,&uPerpLinearUp,OutputPerpLinearUpName+"/F");
  TBranch *uParaLinearUpBr = outTree->Branch(OutputParaLinearUpName,&uParaLinearUp,OutputParaLinearUpName+"/F");
  TBranch *uMagLinearUpBr  = outTree->Branch(OutputRecoilPtLinearUpName,&uMagLinearUp,OutputRecoilPtLinearUpName+"/F");
  TBranch *uPerpLinearDownBr = outTree->Branch(OutputPerpLinearDownName,&uPerpLinearDown,OutputPerpLinearDownName+"/F");
  TBranch *uParaLinearDownBr = outTree->Branch(OutputParaLinearDownName,&uParaLinearDown,OutputParaLinearDownName+"/F");
  TBranch *uMagLinearDownBr  = outTree->Branch(OutputRecoilPtLinearDownName,&uMagLinearDown,OutputRecoilPtLinearDownName+"/F");
  TBranch *uPerpQuadUpBr = outTree->Branch(OutputPerpQuadUpName,&uPerpQuadUp,OutputPerpQuadUpName+"/F");
  TBranch *uParaQuadUpBr = outTree->Branch(OutputParaQuadUpName,&uParaQuadUp,OutputParaQuadUpName+"/F");
  TBranch *uMagQuadUpBr  = outTree->Branch(OutputRecoilPtQuadUpName,&uMagQuadUp,OutputRecoilPtQuadUpName+"/F");
  TBranch *uPerpQuadDownBr = outTree->Branch(OutputPerpQuadDownName,&uPerpQuadDown,OutputPerpQuadDownName+"/F");
  TBranch *uParaQuadDownBr = outTree->Branch(OutputParaQuadDownName,&uParaQuadDown,OutputParaQuadDownName+"/F");
  TBranch *uMagQuadDownBr  = outTree->Branch(OutputRecoilPtQuadDownName,&uMagQuadDown,OutputRecoilPtQuadDownName+"/F");

  int lastPdgId = 0;

  int numEntries = inTree->GetEntriesFast();
  for (int entry = 0; entry != numEntries; ++entry) {
    if (entry % reportFreq == 0)
      std::cout <<  "Processing " << inFileName << " ... " << float(entry)/numEntries * 100 << "%" << std::endl;

    double footPt = 0.;
    double footUp = 0.;
    double footDown = 0.;
    double ZfootPt = 0.;
    double ZfootUp = 0.;
    double ZfootDown = 0.;

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

      ZfootPt   = ZPt + (mmFuncVal - eeFuncVal)/(1 - ZmmFunc->GetParameter(1));
      ZfootUp   = ZPt + (mmFuncVal - eeFuncVal + sqrt(mmUpErr2 + eeDownErr2))/(1 - ZmmFuncUp->GetParameter(1) - 2*ZmmFuncUp->GetParameter(2)*ZPt);
      ZfootDown = ZPt + (mmFuncVal - eeFuncVal - sqrt(mmDownErr2 + eeUpErr2))/(1 - ZmmFuncDown->GetParameter(1) - 2*ZmmFuncDown->GetParameter(2)*ZPt);
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
      rc->ComputeU(genBosPt,uPara,uParaScaleUp,uParaScaleDown,uParaResolutionUp,uParaResolutionDown,
                   uPerp,uPerpResolutionUp,uPerpResolutionDown);
      uMag = TMath::Sqrt(uPerp*uPerp + uPara*uPara);
      uMagScaleUp = TMath::Sqrt(uPerp*uPerp + uParaScaleUp*uParaScaleUp);
      uMagScaleDown = TMath::Sqrt(uPerp*uPerp + uParaScaleDown*uParaScaleDown);
      uMagResolutionUp = TMath::Sqrt(uPerpResolutionUp*uPerpResolutionUp + uParaResolutionUp*uParaResolutionUp);
      uMagResolutionDown = TMath::Sqrt(uPerpResolutionDown*uPerpResolutionDown + uParaResolutionDown*uParaResolutionDown);
      if ((ZfootUp - ZfootDown) > (footUp - footDown)) {
        uParaFootUp = uPara - ZfootUp + ZfootPt;
        uParaFootDown = uPara - ZfootDown + ZfootPt;
      }
      else {
        uParaFootUp = uPara - footUp + footPt;
        uParaFootDown = uPara - footDown + footPt;
      }
      uMagFootUp = TMath::Sqrt(uPerp*uPerp + uParaFootUp*uParaFootUp);
      uMagFootDown = TMath::Sqrt(uPerp*uPerp + uParaFootDown*uParaFootDown);

      uPerpLinearUp = uPerpResolutionUp;
      uParaLinearUp = uParaResolutionUp + uParaFootUp + uParaScaleUp - 2*uPara;
      uMagLinearUp = TMath::Sqrt(uPerpLinearUp*uPerpLinearUp + uParaLinearUp*uParaLinearUp);
      uPerpLinearDown = uPerpResolutionDown;
      uParaLinearDown = uParaResolutionDown + uParaFootDown + uParaScaleDown - 2*uPara;
      uMagLinearDown = TMath::Sqrt(uPerpLinearDown*uPerpLinearDown + uParaLinearDown*uParaLinearDown);

      uPerpQuadUp = uPerpResolutionUp;
      uParaQuadUp = uPara - TMath::Sqrt((uParaResolutionUp - uPara) * (uParaResolutionUp - uPara) +
                                        (uParaFootUp - uPara) * (uParaFootUp - uPara) +
                                        (uParaScaleUp - uPara) * (uParaScaleUp - uPara));
      uMagQuadUp = TMath::Sqrt(uPerpQuadUp*uPerpQuadUp + uParaQuadUp*uParaQuadUp);
      uPerpQuadDown = uPerpResolutionDown;
      uParaQuadDown = uPara + TMath::Sqrt((uParaResolutionDown - uPara) * (uParaResolutionDown - uPara) +
                                        (uParaFootDown - uPara) * (uParaFootDown - uPara) +
                                        (uParaScaleDown - uPara) * (uParaScaleDown - uPara));
      uMagQuadDown = TMath::Sqrt(uPerpQuadDown*uPerpQuadDown + uParaQuadDown*uParaQuadDown);
    }
    else {
      uPerp = 300;
      uPara = 300;
      uMag  = -5;
      uParaScaleUp = 300;
      uMagScaleUp  = -5;
      uParaScaleDown = 300;
      uMagScaleDown  = -5;
      uPerpResolutionUp = 300;
      uParaResolutionUp = 300;
      uMagResolutionUp  = -5;
      uPerpResolutionDown = 300;
      uParaResolutionDown = 300;
      uMagResolutionDown  = -5;
      uParaFootUp = 300;
      uMagFootUp  = -5;
      uParaFootDown = 300;
      uMagFootDown  = -5;
      uPerpLinearUp = 300;
      uParaLinearUp = 300;
      uMagLinearUp  = -5;
      uPerpLinearDown = 300;
      uParaLinearDown = 300;
      uMagLinearDown  = -5;
      uPerpQuadUp = 300;
      uParaQuadUp = 300;
      uMagQuadUp  = -5;
      uPerpQuadDown = 300;
      uParaQuadDown = 300;
      uMagQuadDown  = -5;
    }

    if (inTree == outTree) {
      uPerpBr->Fill();
      uParaBr->Fill();
      uMagBr->Fill();
      uParaScaleUpBr->Fill();
      uMagScaleUpBr->Fill();
      uParaScaleDownBr->Fill();
      uMagScaleDownBr->Fill();
      uPerpResolutionUpBr->Fill();
      uParaResolutionUpBr->Fill();
      uMagResolutionUpBr->Fill();
      uPerpResolutionDownBr->Fill();
      uParaResolutionDownBr->Fill();
      uMagResolutionDownBr->Fill();
      uParaFootUpBr->Fill();
      uMagFootUpBr->Fill();
      uParaFootDownBr->Fill();
      uMagFootDownBr->Fill();
      uPerpLinearUpBr->Fill();
      uParaLinearUpBr->Fill();
      uMagLinearUpBr->Fill();
      uPerpLinearDownBr->Fill();
      uParaLinearDownBr->Fill();
      uMagLinearDownBr->Fill();
      uPerpQuadUpBr->Fill();
      uParaQuadUpBr->Fill();
      uMagQuadUpBr->Fill();
      uPerpQuadDownBr->Fill();
      uParaQuadDownBr->Fill();
      uMagQuadDownBr->Fill();
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
