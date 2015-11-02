#include <iostream>
#include "RecoilCorrector.h"
#include "TVector2.h"

RecoilCorrector::RecoilCorrector() :
  fSingleGaus(false),
  fCurrChannel(kZmm)
{
  rng = new TRandom3();
  inName = "Zmm";
  outName[0] = "Zmm";
  outName[1] = "Zee";
  outName[2] = "Znn";
  outName[3] = "gjets";
  outName[4] = "Wmn";
  outName[5] = "Wen";
  for (unsigned int iR=0; iR!=3; ++iR) {
    for (unsigned int iU=0; iU!=2; ++iU) {
      for (unsigned int iC = 0; iC != 6; ++iC) {
        fmu[iU][iR][iC] = NULL;
        covMu[iU][iR][iC] = NULL;
        fsigma1[iU][iR][iC] = NULL;
        covSigma1[iU][iR][iC] = NULL;
        fsigma2[iU][iR][iC] = NULL;
        covSigma2[iU][iR][iC] = NULL;
        fsigma[iU][iR][iC] = NULL;
        covSigma[iU][iR][iC] = NULL;
        fsigmaSingle[iU][iR][iC] = NULL;
        covSigmaSingle[iU][iR][iC] = NULL;
      }
    }
  } 
}

RecoilCorrector::~RecoilCorrector() { 
  for (unsigned int iR=0; iR!=3; ++iR) {
    for (unsigned int iU=0; iU!=2; ++iU) {
      for (unsigned int iC = 0; iC != 6; ++iC) {
        if (covMu[iU][iR][iC])
          delete covMu[iU][iR][iC];
        if (covSigma1[iU][iR][iC])
          delete covSigma1[iU][iR][iC];
        if (covSigma2[iU][iR][iC])
          delete covSigma2[iU][iR][iC];
        if (covSigma[iU][iR][iC])
          delete covSigma[iU][iR][iC];
        if (covSigmaSingle[iU][iR][iC])
          delete covSigmaSingle[iU][iR][iC];
      }
    } 
  }
  delete rng;
}

void RecoilCorrector::SetFitResult(TF1 *f, TMatrixDSym *cov, RecoilType rType, UType uType, Parameter p, ChannelNum chan) {
  // Store pointer to fit function and recover covariance matrix of fit
  // For now, f cannot be deleted until RecoilCorrector is done running
  // TODO: clone f parameters and store locally, involves assuming poly
  // form of f
  switch (p) {
    case kMu:
      // if (covMu[uType][rType]) {
      //   // delete fmu[uType][rType];
      //   delete covMu[uType][rType][chan];
      // }
      fmu[uType][rType][chan] = f;
      covMu[uType][rType][chan] = (TMatrixDSym*)cov->Clone();
      if (!xxMu)
        xxMu = new double[covMu[uType][rType][chan]->GetNrows()];
      break;
    case kSigma1:
      // if (covSigma1[uType][rType][chan]) {
      //   // delete fsigma1[uType][rType][chan];
      //   delete covSigma1[uType][rType][chan];
      // }
      fsigma1[uType][rType][chan] = f;
      covSigma1[uType][rType][chan] =  (TMatrixDSym*)cov->Clone();
      if (!xxSigma1)
        xxSigma1 = new double[covSigma1[uType][rType][chan]->GetNrows()];
      break;
    case kSigma2:
      // if (covSigma2[uType][rType][chan]) {
      //   // delete fsigma2[uType][rType][chan];
      //   delete covSigma2[uType][rType][chan];
      // }
      fsigma2[uType][rType][chan] = f;
      covSigma2[uType][rType][chan] =  (TMatrixDSym*)cov->Clone();
      if (!xxSigma2)
        xxSigma2 = new double[covSigma2[uType][rType][chan]->GetNrows()];
      break;
    case kSigma:
      // if (covSigma[uType][rType][chan]) {
      //   // delete fsigma[uType][rType][chan];
      //   delete covSigma[uType][rType][chan];
      // }
      fsigma[uType][rType][chan] = f;
      covSigma[uType][rType][chan] =  (TMatrixDSym*)cov->Clone();
      if (!xxSigma)
        xxSigma = new double[covSigma[uType][rType][chan]->GetNrows()];
      break;
    case kSigmaSingle:
      // if (covSigmasingle[uType][rType][chan]) {
      //   // delete fsigmaSingle[uType][rType][chan];
      //   delete covSigmaSingle[uType][rType][chan];
      // }
      fsigmaSingle[uType][rType][chan] = f;
      covSigmaSingle[uType][rType][chan] =  (TMatrixDSym*)cov->Clone();
      if (!xxSigmaSingle)
        xxSigmaSingle = new double[covSigmaSingle[uType][rType][chan]->GetNrows()];
      break;
  }
}

void RecoilCorrector::LoadAllFits(TFile *fIn) {
  TString fitBaseName;
  
  for (unsigned int iC = 0; iC != 6; ++iC) {
    TString recoilNames[3] = {"Data_"+inName,"MC_"+inName,"MC_"+outName[iC]};
    TF1 *f;
    TMatrixDSym *cov;
    // fprintf(stderr,"RecoilCorrector::LoadAllFits: Careful not to close %s until you are done with RecoilCorrector.\n",fIn->GetName());
    for (unsigned int iR=0; iR!=3; ++iR) {
      for (int iU=0; iU!=2; ++iU) {
        fitBaseName = TString::Format("u%i_%s",iU+1,recoilNames[iR].Data());
        // fprintf(stderr,"loading %s\n",fitBaseName.Data()); 
        
        f = (TF1*)fIn->Get("fcn_mu_"+fitBaseName);
        cov = (TMatrixDSym*)fIn->Get("cov_mu_"+fitBaseName);
        SetFitResult(f,cov,(RecoilType)iR,(UType)iU,kMu,(ChannelNum)iC);      
        
        f = (TF1*)fIn->Get("fcn_sig1_"+fitBaseName);
        cov = (TMatrixDSym*)fIn->Get("cov_sig1_"+fitBaseName);
        SetFitResult(f,cov,(RecoilType)iR,(UType)iU,kSigma1,(ChannelNum)iC);      
        
        f = (TF1*)fIn->Get("fcn_sig2_"+fitBaseName);
        cov = (TMatrixDSym*)fIn->Get("cov_sig2_"+fitBaseName);
        SetFitResult(f,cov,(RecoilType)iR,(UType)iU,kSigma2,(ChannelNum)iC);      
        
        f = (TF1*)fIn->Get("fcn_sig3_"+fitBaseName);
        cov = (TMatrixDSym*)fIn->Get("cov_sig3_"+fitBaseName);
        SetFitResult(f,cov,(RecoilType)iR,(UType)iU,kSigma,(ChannelNum)iC);      
        
        f = (TF1*)fIn->Get("fcn_sig_"+fitBaseName);
        cov = (TMatrixDSym*)fIn->Get("cov_sig_"+fitBaseName);
        SetFitResult(f,cov,(RecoilType)iR,(UType)iU,kSigmaSingle,(ChannelNum)iC);      
      }
    } // loop over u1 u2
  } // loop over recoil types
}

double RecoilCorrector::GetError(double x,RecoilType r,UType u,Parameter p,ChannelNum c) const {
  TMatrixDSym *cov = 0;
  double *xx = 0;
  switch (p) {
    case kMu:
      cov = covMu[u][r][c];
      xx = xxMu;
      break;
    case kSigma1:
      cov = covSigma1[u][r][c];
      xx = xxSigma1;
      break;
    case kSigma2:
      cov = covSigma2[u][r][c];
      xx = xxSigma2;
      break;
    case kSigma:
      cov = covSigma[u][r][c];
      xx = xxSigma;
      break;
    case kSigmaSingle:
      cov = covSigmaSingle[u][r][c];
      xx = xxSigmaSingle;
      break;
  }

  unsigned int dim = cov->GetNrows(); 
  
  for (unsigned int iR=0; iR!=dim; ++dim) 
    xx[iR] = TMath::Power(x,(Int_t)iR);

  double error=0;
  for (unsigned int iR=0; iR!=dim; ++iR) {
    error += xx[iR] * xx[iR] * (*cov)(iR,iR);
    for (unsigned int iC=0; iC!=iR; ++iC) {
      error += 2* xx[iC] * (*cov)(iC,iR) * xx[iR];
    }
  }

  return error;
}

void RecoilCorrector::ComputeU(float genpt, float &u1, float &u2, float nsigma/*=0*/) const {
  
  // first compute u1
  //double mu     = (fmu[kU1][kDataIn][0]->Eval(genpt))     * (fmu[kU1][kMCOut][fCurrChannel]->Eval(genpt))     / (fmu[kU1][kMCIn][0]->Eval(genpt));
  double mu     = (fmu[kU1][kDataIn][0]->Eval(genpt)-genpt)     * (fmu[kU1][kMCOut][fCurrChannel]->Eval(genpt)-genpt)     / (fmu[kU1][kMCIn][0]->Eval(genpt)-genpt);
  double sigma1 = fsigma1[kU1][kDataIn][0]->Eval(genpt) * fsigma1[kU1][kMCOut][fCurrChannel]->Eval(genpt) / fsigma1[kU1][kMCIn][0]->Eval(genpt);
  double sigma2 = fsigma2[kU1][kDataIn][0]->Eval(genpt) * fsigma2[kU1][kMCOut][fCurrChannel]->Eval(genpt) / fsigma2[kU1][kMCIn][0]->Eval(genpt);
  double sigma  = fsigma[kU1][kDataIn][0]->Eval(genpt)  * fsigma[kU1][kMCOut][fCurrChannel]->Eval(genpt)  / fsigma[kU1][kMCIn][0]->Eval(genpt);
  double sigmaSingle  = fsigmaSingle[kU1][kDataIn][0]->Eval(genpt)  * fsigmaSingle[kU1][kMCOut][fCurrChannel]->Eval(genpt)  / fsigmaSingle[kU1][kMCIn][0]->Eval(genpt);

  // a la error propogation used in w/z analaysis
  // TODO: improve to treat parameters independently
  // currently making conservative assumption of maximal correlation
  if (nsigma != 0.) {
    mu     += nsigma * GetError(genpt,kMCOut,kU1,kMu,fCurrChannel)     * (fmu[kU1][kDataIn][0]->Eval(genpt)-genpt) / (fmu[kU1][kMCIn][0]->Eval(genpt) - genpt);
    sigma1 += nsigma * GetError(genpt,kMCOut,kU1,kSigma1,fCurrChannel) * fsigma1[kU1][kDataIn][0]->Eval(genpt) / fsigma1[kU1][kMCIn][0]->Eval(genpt); 
    sigma2 += nsigma * GetError(genpt,kMCOut,kU1,kSigma2,fCurrChannel) * fsigma2[kU1][kDataIn][0]->Eval(genpt) / fsigma2[kU1][kMCIn][0]->Eval(genpt);
    sigma  += nsigma * GetError(genpt,kMCOut,kU1,kSigma,fCurrChannel)  * fsigma[kU1][kDataIn][0]->Eval(genpt)  / fsigma[kU1][kMCIn][0]->Eval(genpt);
    sigmaSingle  += nsigma * GetError(genpt,kMCOut,kU1,kSigmaSingle,fCurrChannel)  * fsigmaSingle[kU1][kDataIn][0]->Eval(genpt)  / fsigmaSingle[kU1][kMCIn][0]->Eval(genpt);
  }

  double frac = (sigma-sigma2)/(sigma1-sigma2);

  if (fSingleGaus)
    u1 = rng->Gaus(mu,sigmaSingle);
  else
    u1 = (((rng->Uniform(0,1)<frac) ? rng->Gaus(mu,sigma1) : rng->Gaus(mu,sigma2)));

  // now compute u2
  mu     = fmu[kU2][kDataIn][0]->Eval(genpt)     * fmu[kU2][kMCOut][fCurrChannel]->Eval(genpt)     / fmu[kU2][kMCIn][0]->Eval(genpt);
  sigma1 = fsigma1[kU2][kDataIn][0]->Eval(genpt) * fsigma1[kU2][kMCOut][fCurrChannel]->Eval(genpt) / fsigma1[kU2][kMCIn][0]->Eval(genpt);
  sigma2 = fsigma2[kU2][kDataIn][0]->Eval(genpt) * fsigma2[kU2][kMCOut][fCurrChannel]->Eval(genpt) / fsigma2[kU2][kMCIn][0]->Eval(genpt);
  sigma  = fsigma[kU2][kDataIn][0]->Eval(genpt)  * fsigma[kU2][kMCOut][fCurrChannel]->Eval(genpt)  / fsigma[kU2][kMCIn][0]->Eval(genpt);
  sigmaSingle  = fsigmaSingle[kU2][kDataIn][0]->Eval(genpt)  * fsigmaSingle[kU2][kMCOut][fCurrChannel]->Eval(genpt)  / fsigmaSingle[kU2][kMCIn][0]->Eval(genpt);
  
  if (nsigma != 0.) {
    mu     += nsigma * GetError(genpt,kMCOut,kU2,kMu,fCurrChannel)     * fmu[kU2][kDataIn][0]->Eval(genpt)     / fmu[kU2][kMCIn][0]->Eval(genpt);
    sigma1 += nsigma * GetError(genpt,kMCOut,kU2,kSigma1,fCurrChannel) * fsigma1[kU2][kDataIn][0]->Eval(genpt) / fsigma1[kU2][kMCIn][0]->Eval(genpt); 
    sigma2 += nsigma * GetError(genpt,kMCOut,kU2,kSigma2,fCurrChannel) * fsigma2[kU2][kDataIn][0]->Eval(genpt) / fsigma2[kU2][kMCIn][0]->Eval(genpt);
    sigma  += nsigma * GetError(genpt,kMCOut,kU2,kSigma,fCurrChannel)  * fsigma[kU2][kDataIn][0]->Eval(genpt)  / fsigma[kU2][kMCIn][0]->Eval(genpt);
    sigmaSingle  += nsigma * GetError(genpt,kMCOut,kU2,kSigmaSingle,fCurrChannel)  * fsigmaSingle[kU2][kDataIn][0]->Eval(genpt)  / fsigmaSingle[kU2][kMCIn][0]->Eval(genpt);
  }

  frac = (sigma-sigma2)/(sigma1-sigma2);

  if (fSingleGaus)
    u2 = rng->Gaus(mu,sigmaSingle);
  else
    u2 = (rng->Uniform(0,1)<frac) ? rng->Gaus(mu,sigma1) : rng->Gaus(mu,sigma2);
}

void RecoilCorrector::CorrectMET(float genpt,float genphi,float leppt,float lepphi,float& met,float& metphi, float nsigma, float u1, float u2) const {
  if (u1==-999||u2==-999)
    ComputeU(genpt,u1,u2,nsigma);
  TVector2 vLep(leppt*TMath::Cos(lepphi),leppt*TMath::Sin(lepphi));
  TVector2 vU(u1*TMath::Cos(genphi)-u2*TMath::Sin(genphi),u1*TMath::Sin(genphi)+u2*TMath::Cos(genphi));
  TVector2 vMissingEnergy = -1*(vLep+vU);
  met = vMissingEnergy.Mod();
  metphi = vMissingEnergy.Phi();
}
