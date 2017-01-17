#ifndef KMEANSWEIGHTS_H
#define KMEANSWEIGHTS_H

class kMeansWeights
{
 public: 
  kMeansWeights(Int_t k) ;
  ~kMeansWeights() {}
  void fit();
  void ScatterPlot();
  void makeHistos();
  void MCScatterPlots();
 protected:
  Int_t      fK;     // numer of clusters
  Double_t*  fMx;    // centroid position X
  Double_t*  fMy;    // centroid position Y
  vector<Double_t>  fX;     // data component X
  vector<Double_t>  fY;     // data component Y
  vector<Double_t>  fW;     // data weight

  Int_t*     fTgt;   // target
  Int_t      fN;     // number of entries
  vector<Double_t>* fDataX; // data in clusters
  vector<Double_t>* fDataY; // data in clusters

  // Monte Carlo
  vector<Double_t>  fTTbarX;     // data component X
  vector<Double_t>  fTTbarY;     // data component Y
  vector<Double_t>  fTTbarW;     // data weight
  vector<Double_t>  fTTHX;     // data component X
  vector<Double_t>  fTTHY;     // data component Y
  vector<Double_t>  fTTHW;     // data weight
  vector<Double_t>  fTTWX;     // data component X
  vector<Double_t>  fTTWY;     // data component Y
  vector<Double_t>  fTTWW;     // data weight

  void     Init();
  Double_t d(Double_t, Double_t, Double_t, Double_t);
  void     MakeAssignment();
  Int_t    CalculateCentroids();
  void     readFromFiles();
  Int_t    GetCluster(Double_t, Double_t);
  Int_t    classicalBinning(Double_t, Double_t);
};


#endif
