#ifndef MAKESIMPLECARD_H
#define MAKESIMPLECARD_H

#include <vector>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TString.h>
#include <TFile.h>

using namespace std;

class MakeSimpleCard
{
 public:
  MakeSimpleCard(TH1* sig, vector<TH1*> bkg, TString cardName="datacard.txt", double lumi=37000., bool debug=false);
  ~MakeSimpleCard() {};
  
  void doCard();
  
 protected:

  void Renormalize();
  void FillFirstBlock();
  void FillObserved();
  void FillRates();
  void WriteShapes();
  void DoStatVariation(TH1* shape, TString basename);
  void WriteCard();

  bool debug_;
  double lumi_;

  TH1* sig_;
  vector<TH1*> bkg_;
  TH1* data_;
  
  vector<TH1*> statVariations_;

  TString cardName_;

  TFile* shapesFile_;

  double systSig_;
  double systBkg_;
  
  size_t nBkg_;
  
  ofstream card_;

};

#endif
