#ifndef RECURSIVECLUSTERING_H
#define RECURSIVECLUSTERING_H
#include <utility>

#include"TString.h"
#include"TFormula.h"

class Point
{
 public:
  Double_t fX;
  Double_t fY;
  Double_t fW;
  Point(){}
  Point(Double_t x, Double_t y, Double_t w){ fX = x; fY = y; fW = w;}

};


class Cluster
{
 public:
  Cluster(){}
  Cluster( vector<Point> , vector<Point>, vector<Point>, Int_t, TString );
  ~Cluster() { } 
  void recluster();
  Int_t    FindUnclusterizableCluster( Point );
  TString  FindSubClusterName(Point, Int_t);
  TString  GetName(){ return fName; }

 private:
  vector<Point> fTTbar;
  vector<Point> fTTH  ;
  vector<Point> fTTW  ;
  vector<Point> fData ;

  vector<Point> fCentroids;

  vector<Int_t> fTgt;
  vector<Cluster> SubClusters;
  TString fName; 
  TRandom3* r_;

  bool fIsClusterizable;
  Int_t fIndex; // global label for every non clusterizable cluster

  void     MakeAssignment();
  Int_t    CalculateCentroids();
  UInt_t    fK;
  UInt_t    FindSubCluster( Point, bool verbose=false );
  Double_t d(Double_t, Double_t, Double_t, Double_t);

  double FOMforClusterCut( const double*);
  TFormula fFormula;


};

// https://english.stackexchange.com/questions/205815/is-targetted-a-standard-british-english-spelling
class TargettedClustering
{
 public: 
  TargettedClustering(Int_t k, Int_t nLep=2, UInt_t seed=0) ;
  ~TargettedClustering() {}
  void makeHistos();
  void VoronoiPlot();
  void VoronoiPlot2(Int_t);
  void     Test();
  Double_t SignificanceAtLevel(Int_t);
  void SignificancesEachLevel();

 protected:
  // Monte Carlo
  vector<Point>  fTTbar;  
  vector<Point>  fTTH;    
  vector<Point>  fTTW;    
  vector<Point> fTTbarMC;
  vector<Point> fTTHMC;
  vector<Point> fTTWMC;
  vector<Point> fCentroids;
  vector<double> fSignificances;
  Int_t fK;
  Int_t nLep_; // Choose between 2lss and 3l input
  UInt_t seed_;

  TRandom3* r_; // Randomizer

  Cluster mainCluster;
  Double_t d(Double_t, Double_t, Double_t, Double_t);
  void     readFromFiles();
  void     StartTheThing();
  void     StoreToFile();
};


#endif
