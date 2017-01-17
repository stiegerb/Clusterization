#include "RecursiveClustering.h"
#include <iostream>

using namespace std;

void Cluster::ReadEvents(TString& fNameTTH, TString& fNameBkg)
{
  // Fixme: make fNameBkg a vector.
  fOut_ = TFile::Open("events.root", "RECREATE");

  TTree *tTTH_ = new TTree("tth", "Data from tth signal");
  TTree *tBkg_ = new TTree("bkg", "Data from bkg background");

  Long64_t nTTH(tTTH_->ReadFile(fNameTTH,"bdt1:bdt2:weight"));
  Long64_t nBkg(tBkg_->ReadFile(fNameBkg,"bdg1:bdt2:weight"));
  
  cout << "[Cluster::ReadEvents] Events read:" << endl;
  cout << "[Cluster::ReadEvents] \t tth: " << nTTH << endl;
  cout << "[Cluster::ReadEvents] \t bkg: " << nBkg << endl;
  
  tTTH_->Write();
  tBkg_->Write();

}

Cluster::Cluster(double* center, TString fNameTTH, TString fNameBkg)
{
  // Fixme: make fNameBkg a vector.
  ReadEvents(fNameTTH,fNameBkg);
  Init();

}


void Cluster::Init()
{
}

void Cluster::Recluster()
{
}

void Cluster::FindSubCluster(double* point)
{
}


RecursiveClustering::RecursiveClustering()
{


}


