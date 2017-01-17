#ifndef RECURSIVECLUSTERING_H
#define RECURSIVECLUSTERING_H

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

class Cluster
{
 public:
  Cluster(double* center, TString fNameTTH="data/tth.txt", TString fNameBkg="data/ttbar.txt");
  ~Cluster() {};

 protected:
  void ReadEvents(TString& fNameTTH, TString& fNameBkg);
  void Init();
  void Recluster();
  void FindSubCluster(double* point);

  double* center_;
  TFile* fOut_; // Not sure it's needed. Perhaps to avoid memory-resident trees
  TTree* tTTH_;
  TTree* tBkg_;
};


class RecursiveClustering
{
 public:
  RecursiveClustering();
  ~RecursiveClustering() {};
  void GetEventsFromFile();
  void clusterize();
  void makeHistos();
  void makeScatterPlot();

  TFile* fOut_; // Not sure it's needed. Perhaps to avoid memory-resident trees

  TTree*
    trainTTbar_,
    trainTTV_,
    trainTTH,
    testTTbar,
    testTTV,
    testTTH;
  
};

#endif
