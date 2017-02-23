#ifndef LikeliHOOD_H
#define LikeliHOOD_H
#include <utility>

#include"TString.h"
#include"TGraph.h"
#include"TH2.h"

class Point
{
 public:
  Double_t fX;
  Double_t fY;
  Double_t fW;
  Point(){}
  Point(Double_t x, Double_t y, Double_t w){ fX = x; fY = y; fW = w;}

};

class LikeliHOOD
{
 public:
  LikeliHOOD( Int_t nLep=2, TString fileType="txt", Int_t trial=0) ;
  ~LikeliHOOD() {}
  // void VoronoiPlot();
  //  void Test();
  void SortedTest();
  void Recluster();
  void Test();
  Double_t GetFOM(Int_t, Int_t);
  //  void ReCleanClusters();
  void LoadClusterFromFile();
  void     StoreToFile();
  void     VoronoiPlot();

 protected:
  // Monte Carlo
  vector<Point>  fTTbar;
  vector<Point>  fTTH;
  vector<Point>  fTTW;
  vector<Point> fTTbarMC;
  vector<Point> fTTHMC;
  vector<Point> fTTWMC;
  vector<TGraph*> contours;
  void MakeLikeliHood();
  // vector<vector<Int_t>> clusters; // each array contains an array of all the miniclusters that it is made of
  Int_t nLep_; // Choose between 2lss and 3l input

  Int_t trial_;
  vector<pair<Double_t,Int_t>> SoverB;
  void readFromTxTFiles();
  void readFromRootFiles();
  Double_t GetCluster(Point);
  TString fileType;
  //  Int_t   NearestClusters(Int_t);
  void     Init();
  void     ReadFromFiles();
  void     StartTheThing();
  Int_t    SortedThing(Int_t);
  void     MakeFineBinning();
  void     ReMakeTarget();
  void     GetTheContours();
  void     GETCUM();
  Double_t GetLikeLiHood(Point);

  TH2F*    hBkg;
  TH2F*    hSig;
  TH2F*    hRtio;
  TH2F*    hTargetBinning;
  TH1F*    hFineTTbar;
  TH1F*    hFineTTW;
  TH1F*    hFineTTH;
  TH1F*    hAuxHisto;
};


#endif
