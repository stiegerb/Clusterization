#include "RecursiveClustering.h"
#include <iostream>

using namespace std;

void Cluster::ReadEvents(TString& fNameTTH, TString& fNameBkg)
{
  // Fixme: make fNameBkg a vector.
  fOut_ = TFile::Open("events.root", "RECREATE");

  tTTH_ = new TTree("tth", "Data from tth signal");
  tBkg_ = new TTree("bkg", "Data from bkg background");

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

void RecursiveClustering::GetEventsFromFile()
{

  // Fixme: make fNameBkg a vector.
  fOut_ = TFile::Open("events.root", "RECREATE");

  TTree* tTTb = new TTree("ttbar","Data from ttbar");
  TTree* tTTV = new TTree("ttv","Data from ttV");
  TTree* tTTH = new TTree("tth","Data from ttH");

  Long64_t nTTb(tTTb->ReadFile("data/ttbar.txt","bdt1:bdt2:weight"));
  Long64_t nTTV(tTTV->ReadFile("data/ttv.txt"  ,"bdg1:bdt2:weight"));
  Long64_t nTTH(tTTH->ReadFile("data/tth.txt"  ,"bdt1:bdt2:weight"));

  cout << "[RecursiveClustering::GetEventsFromFile] Events read:" << endl;
  cout << "[RecursiveClustering::GetEventsFromFile] \t tth:   " << nTTH << endl;
  cout << "[RecursiveClustering::GetEventsFromFile] \t ttv:   " << nTTV << endl;
  cout << "[RecursiveClustering::GetEventsFromFile] \t ttbar: " << nTTb << endl;

  // Now split into training and testing with clonetree

  tTTH->Write();
  tBkg->Write();

}
