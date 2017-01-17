#ifndef RECURSIVECLUSTERING_H
#define RECURSIVECLUSTERING_H
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
  Cluster( vector<Point> , vector<Point>, vector<Point>, Int_t  );
  ~Cluster() { } 
  void recluster();
  Int_t    FindUnclusterizableCluster( Point );
  
 private:
  vector<Point> fTTbar;
  vector<Point> fTTH  ;
  vector<Point> fTTW  ;
  vector<Point> fData;
  
  vector<Point> fCentroids;
  vector<Int_t> fTgt;
  vector<Cluster> SubClusters;

  bool fIsClusterizable;
  Int_t fIndex; // global label for every non clusterizable cluster

  void     MakeAssignment();
  Int_t    CalculateCentroids();
  Int_t    fK;
  Int_t    FindSubCluster( Point );
  Double_t d(Double_t, Double_t, Double_t, Double_t);
};

class RecursiveClustering
{
 public: 
  RecursiveClustering(Int_t k) ;
  ~RecursiveClustering() {}
  void fit();
  //void ScatterPlot();
  void makeHistos();
  void     Test();

 protected:
  // Monte Carlo
  vector<Point>  fTTbar;  
  vector<Point>  fTTH;    
  vector<Point>  fTTW;    
  vector<Point> fTTbarMC;
  vector<Point> fTTHMC;
  vector<Point> fTTWMC;
  Int_t fK;

  Cluster mainCluster;
  void     Init();
  Double_t d(Double_t, Double_t, Double_t, Double_t);
  void     readFromFiles();
  void     StartTheThing();
};


#endif
