#ifndef KMEANSWEIGHTS_H
#define KMEANSWEIGHTS_H

// Forward declarations
#include <TString.h> 

class kMeansWeights
{
 public: 
  kMeansWeights(Int_t k) ;
  ~kMeansWeights() {}
  void fit();
  void ScatterPlot();

 protected:
  Int_t      fK;     // number of clusters
  TString     fKname; // number of clusters (string)
  Double_t*  fMx;    // centroid position X
  Double_t*  fMy;    // centroid position Y
  vector<Double_t>  fX;     // data component X
  vector<Double_t>  fY;     // data component Y
  vector<Double_t>  fW;     // data weight
  Int_t*     fTgt;   // target
  Int_t      fN;     // number of entries
  vector<Double_t>* fDataX; // data in clusters
  vector<Double_t>* fDataY; // data in clusters

  void     Init();
  Double_t d(Double_t, Double_t, Double_t, Double_t);
  void     MakeAssignment();
  Int_t    CalculateCentroids();
  void     readFromFiles();

};


#endif
