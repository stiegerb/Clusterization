#ifndef BOTTOMUPCLUSTERING_H
#define BOTTOMUPCLUSTERING_H
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

class BottomUpClustering
{
 public:
  BottomUpClustering( Int_t nLep=2, TString fileType="txt", Int_t trial=0) ;
  ~BottomUpClustering() {}
  void VoronoiPlot();
  //  void Test();
  void SortedTest();
  void Recluster();

 protected:
  // Monte Carlo
  vector<Point>  fTTbar;
  vector<Point>  fTTH;
  vector<Point>  fTTW;
  vector<Point> fTTbarMC;
  vector<Point> fTTHMC;
  vector<Point> fTTWMC;
  vector<TGraph*> contours;

  vector<vector<Int_t>> clusters; // each array contains an array of all the miniclusters that it is made of
  Int_t nLep_; // Choose between 2lss and 3l input

  Int_t trial_;
  vector<pair<Double_t,Int_t>> SoverB;
  void readFromTxTFiles();
  void readFromRootFiles();

  TString fileType;

  void     Init();
  void     ReadFromFiles();
  void     StartTheThing();
  void     StoreToFile();
  Int_t    SortedThing(Int_t);
  bool     AreClustersTogether(int,int);
  void     MakeFineBinning();
  void     ReMakeTarget();
  void     GetEventsInCluster(Int_t, Double_t&, Double_t&, Double_t&, Double_t&, Double_t&, Double_t&);
  Double_t GetFOM(Int_t, Int_t);
  void     GetTheContours();

  TH2F*    hFineBinning;
  TH2F*    hTargetBinning;
};


#endif
