#ifndef RECOILCORRECTOR
#define RECOILCORRECTOR 1

#include <vector>
#include "TMath.h"
#include <assert.h>
#include "TF1.h"
#include "TFile.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TString.h"

/*
Simple class that corrects MET given a pt-binned double Gaussian fit 
to U1 and U2 in 3 samples: input(data), input(MC), and the sample that 
needs to be corrected (MC). Can either compute U1,U2 given boson pT, or 
can correct MET given boson and lepton transverse vectors.
Author(s): S.Narayanan
*/

class RecoilCorrector
{
public:
  RecoilCorrector();
  ~RecoilCorrector();

  enum RecoilType {
    kDataIn,
    kMCIn,
    kMCOut
  };

  enum CorrectionType {
    kType1,
    kType2, // not implemented yet
    kTypeAll
  };

  enum UType {
    kU1,
    kU2
  };

  enum ChannelNum {
    kZmm = 0,
    kZee,
    kZnn,
    kGJets,
    kWmn,
    kWen
  };

  enum Parameter {
    kMu,
    kSigma1,
    kSigma2,
    kSigma,
    kSigmaSingle
  };

  // initialize
  void SetCorrectionType(CorrectionType k)  { fCorrectionType = k; }
  void SetInputName(const char *s)          { inName = s;          }
  void SetOutputName(const char *s)         { outName = s;         } 
  void SetSeed(unsigned int u)              { rng->SetSeed(u);     }

  void SetOutput(ChannelNum c)              { fCurrChannel = c;    }

  // read inputs
  void SetFitResult(TF1*,TMatrixDSym*,RecoilType,UType,Parameter);
  void LoadAllFits(TFile*);

  // computations
  double GetError(double,RecoilType,UType,Parameter) const;
  void ComputeU(float genpt,float& u1,float& u2, float nsigma=0) const;
  void CorrectMET(float genpt,float genphi,float leppt,float lepphi,float& met,float& metphi, float nsigma=0, float u1=-999, float u2=-999) const;

  void SetSingleGaus (bool single) { fSingleGaus = single; }

protected:
  // unsigned int nMuParams;
  // unsigned int nSigma1Params;
  // unsigned int nSigma2Params;
  // unsigned int nSigmaParams;

  Bool_t fSingleGaus;

  CorrectionType fCorrectionType;

  // in principle could get the TF1s from the fit results
  // but this way, we can save everything in a file
  TF1 *fmu[2][3][6]; // U[1,2] mean as a function of pT for Z(data), Z(mc), thing we are trying to correct
  TF1 *fsigma1[2][3][6]; // sigma1
  TF1 *fsigma2[2][3][6]; // sigma2
  TF1 *fsigma[2][3][6]; // sigma
  TF1 *fsigmaSingle[2][3][6];

  // directly store covariance matrices instead of fit results
  TMatrixDSym *covMu[2][3][6];
  TMatrixDSym *covSigma1[2][3][6];
  TMatrixDSym *covSigma2[2][3][6];
  TMatrixDSym *covSigma[2][3][6];
  TMatrixDSym *covSigmaSingle[2][3][6];

  TRandom3 *rng;

  // arrays used for computing errors - allocate them once
  double *xxMu{0};
  double *xxSigma1{0};
  double *xxSigma2{0};
  double *xxSigma{0};
  double *xxSigmaSingle{0};

  TString inName;
  TString outName[6];

  ChannelNum fCurrChannel;

};

#endif
