#include "RecursiveClustering.h"
#include "MakeSimpleCard.h"
#ifdef SIGNIFICANCE_H
#include "Significance.h"
#endif
#include <iostream>
#include <fstream>
#include <vector>
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"
#include <sstream>
#include <TString.h> 
#include "TStyle.h"

using namespace std;

bool SortX( Point i, Point j)
{
  return ( i.fX < j.fX );
}
bool SortY( Point i, Point j)
{
  return ( i.fY < j.fY );
}

std::ostream &operator<<(std::ostream &os, Point const &point) { 
  return os << "(" << point.fX << ", " << point.fY << ")";
}
Int_t gIndex = 0;

Cluster::Cluster( vector<Point> TTbar, vector<Point> TTH  , vector<Point> TTW, Int_t k , TString name, Point centroid)
{
  fName = name;
  cout << "Initialising cluster with name " << name << endl;
  fTTbar  = TTbar;
  fTTH    = TTH  ;	
  fTTW    = TTW  ;	

  fData = fTTbar; 
  fData.insert( fData.end(), fTTW.begin(),fTTW.end()); 
  fData.insert( fData.end(), fTTH.begin(),fTTH.end());

  fIsClusterizable = true;
  fK = k;
  fCentroid = centroid;

}

void Cluster::MakeAssignment()
{
  fTgt.clear();
  // Assigns each data point to its cluster :)
  for (size_t n = 0; n < fData.size(); ++n){
    fTgt.push_back( FindSubCluster( fData[n] ));
  }
}

UInt_t Cluster::FindSubCluster( Point point , bool verbose)
{
  if (verbose ) cout << "Finding subcluster" << endl;
  Int_t cluster = -1;
  Double_t dist = 10000;
  for (size_t k = 0; k < fCentroids.size(); ++k){
    if (verbose ) cout << "Checking centroid " << fCentroids[k] << endl;
    Double_t dst = d( point.fX, point.fY, fCentroids[k].fX, fCentroids[k].fY);
    if (dst < dist){
      cluster = k;
      dist = dst;
    }
  }
  if (cluster < 0){ cout << "[E]: Nearest cluster not found" << endl;}
  return cluster; 
}

Int_t Cluster::FindUnclusterizableCluster( Point point )
{
  // cout << "################################################" << endl;
  // cout << "################################################" << endl;
  // cout << "#########FindUnclusterizableCluster#############" << endl;
  // cout << "################################################" << endl;
  // cout << "################################################" << endl;
  // cout << "#### Cluster " << fName << "################" << endl;
  // cout << point << endl;


  if (fIsClusterizable){
    //    cout << "Checking point in subclusters..." << endl;
    return SubClusters[FindSubCluster(point,false)].FindUnclusterizableCluster(point);
  }
  else{
    return fIndex;
  }

}


TString Cluster::FindSubClusterName( Point point, Int_t level )
{
  if ( fName.Sizeof() <= level and fIsClusterizable){
    return SubClusters[FindSubCluster(point,false)].FindSubClusterName(point,level);
  }
  else{
    return fName;
  }

}

Int_t Cluster::CalculateCentroids()
{
  Int_t theyChange = 1;
  vector<Double_t> fMWght;
  vector<Point> fOldCentroids = fCentroids;

  for (size_t k = 0; k < fCentroids.size(); ++k){
    fCentroids[k].fX = 0;
    fCentroids[k].fY = 0;
    fCentroids[k].fW = -1000000; // weights should never be used for centroids
    fMWght.push_back(0);
  }

  for (size_t n = 0; n < fData.size(); ++n){
    fCentroids[fTgt[n]].fX    +=  fData[n].fW * fData[n].fX;
    fCentroids[fTgt[n]].fY    +=  fData[n].fW * fData[n].fY;
    fMWght[fTgt[n]] +=  fData[n].fW;
  }

  Double_t di = 0;
  
  for (size_t k = 0; k < fCentroids.size(); ++k){
    fCentroids[k].fX /= fMWght[k];
    fCentroids[k].fY /= fMWght[k];
    di += d(fCentroids[k].fX,fCentroids[k].fY, fOldCentroids[k].fX, fOldCentroids[k].fY);
  }
  if (di < 1.0e-8) theyChange =  -1;
  // cout << "[I]: New centroids are (" << fMx[0] << "," << fMy[0]
  //      << ") and (" << fMx[1] << "," << fMy[1] << ")" << endl;


  return theyChange;

}


vector<Point> Cluster::recluster()
{
  TRandom3* r = new TRandom3();
  Double_t maxX = (*max_element(fData.begin(), fData.end(),SortX)).fX;
  Double_t minX = (*min_element(fData.begin(), fData.end(),SortX)).fX;
  Double_t maxY = (*max_element(fData.begin(), fData.end(),SortY)).fY;
  Double_t minY = (*min_element(fData.begin(), fData.end(),SortY)).fY;

  // Random inizialitation of centroids
  for (unsigned int k = 0; k < fK; ++k){
    Point centroid( r->Uniform(minX, maxX), r->Uniform(minY, maxY), -1);
    fCentroids.push_back(centroid);
  }
  // Cluster assignation
  Int_t maxIt = 999999;
  Int_t it = 0;

  while (it < maxIt){
    it++;
    MakeAssignment();
    if (CalculateCentroids() < 0) break;

  }
  cout << "Centroids are " << "(" << fCentroids[0].fX << "," << fCentroids[0].fY 
       << ")" << endl;
  if ( maxIt == it) cout << "[W]: Maximum number of iterations reached!! ";
  
  
  vector<Point> subCentroidList; subCentroidList.clear();
  for (unsigned int k = 0; k < fK; ++k){
    vector<Point> TTbar; TTbar.clear();
    vector<Point> TTW;   TTW.clear();
    vector<Point> TTH;   TTH.clear();
    for (size_t n = 0; n < fTTH.size(); ++n){
      if ( FindSubCluster( fTTH[n] ) == k) TTH.push_back( fTTH[n]);
    }
    for (size_t n = 0; n < fTTW.size(); ++n){
      if ( FindSubCluster( fTTW[n] ) == k) TTW.push_back( fTTW[n]);
    }
    for (size_t n = 0; n < fTTbar.size(); ++n){
      if ( FindSubCluster( fTTbar[n] ) == k) TTbar.push_back( fTTbar[n]);
    }
    if (TTH.size() == 0) fIsClusterizable = false;
    else if (TTbar.size() == 0) fIsClusterizable = false;
    else if (TTW  .size() == 0) fIsClusterizable = false;
    else{
      if (37000*TTH.size()*TTH[0].fW < 5.)
	fIsClusterizable = false;
    }
    Double_t tth = 0;
    Double_t ttw = 0;
    Double_t ttbar = 0;
    if (TTH.size() > 0)   tth   = 36400 * TTH.size()  * TTH[0]  .fW;
    if (TTW.size() > 0)   ttw   = 36400 * TTW.size()  * TTW[0]  .fW;
    if (TTbar.size() > 0) ttbar = 36400 * TTbar.size()  * TTbar[0]  .fW;
    cout << "Expected events in the subcluster " << tth
	 << " " << ttw << " " << ttbar << endl;

    if (!fIsClusterizable){
      cout << "Not clusterizable anymore " << endl;
      cout << "This will be subcluster " << gIndex <<  endl;
      cout << "The centroid is " << fCentroid << endl;
      cout << "The name is " << fName << endl;
      SubClusters.clear();
      fIndex = gIndex;
      gIndex++;
      subCentroidList.clear();
      subCentroidList.push_back( fCentroid );
      break;
    }
    else{
      cout << "Apparently subcluster " << k << " from " << fName 
	   << " is huge!!! (thats what she said), so theres not showstopper not to keep clustering" << endl;
      Cluster subCluster = Cluster( TTbar,  TTH  ,  TTW, fK, fName + Form("%u",k), fCentroids[k]);
      SubClusters.push_back(subCluster);
    }
  }
  for (size_t k = 0; k < SubClusters.size(); ++k){
      cout << "Cluster " << fName << " is huge!!! (thats what she said)" << endl;
      vector<Point> listOfSubcentroids = SubClusters[k].recluster();
      subCentroidList.insert( subCentroidList.end(), listOfSubcentroids.begin(), listOfSubcentroids.end());
  }


  delete r;
  cout << "Name and size of centroid list " << fName << " " << subCentroidList.size() << endl;
  return subCentroidList;

}


RecursiveClustering::RecursiveClustering(Int_t k)
{
  fK = k;
  readFromFiles();
  StartTheThing();
  
  cout << "Produced " << gIndex << " clusters" << endl;
  Point point(0.353811, 0.456623,-1.);
  cout << "Point is " << mainCluster.FindUnclusterizableCluster(point) <<endl;
}

void RecursiveClustering::StartTheThing()
{
  mainCluster = Cluster(fTTbar, fTTH, fTTW, fK, "A", Point(9999,99999,-1));
  fCentroids = mainCluster.recluster();
  StoreToFile();
}

void RecursiveClustering::readFromFiles()
{
  
  fTTbar.clear(); fTTH.clear(); fTTW.clear();
  ifstream f; f.open("data/TTSingleLeptonMarco.txt");
  Int_t count = 0;
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTbar.push_back(point);
    else
      fTTbarMC.push_back(point);
    count++;
  }
  f.close();
  
  f.open("data/ttH.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTH.push_back(point);
    else
      fTTHMC.push_back(point);
    count++;
  }
  f.close();

  f.open("data/ttw.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTW.push_back(point);
    else
      fTTWMC.push_back(point);
    count++;
  }
  f.close();

}

// // Double_t RecursiveClustering::getMax( Double_t* data, Double_t coeff)
// // {
// //   Double_t max = -9999;
// //   for (int n = 0; n < fN; ++n){
// //     if ( coeff*data[n] > max)
// //       max = data[n];
// //   }
// //   if ( max == -9999) cout << "[E]: Maximum hasnt been calculated" << endl;
// //   return max;
// // }

Double_t Cluster::d(Double_t x, Double_t y, Double_t m_x, Double_t m_y) 
{
  return TMath::Sqrt( (x-m_x) * (x-m_x) + (y-m_y) * (y-m_y) );
  //return TMath::Abs(x-m_x) + TMath::Abs(y-m_y);
}

void RecursiveClustering::Test()
{
  TH1F* hTTbar = new TH1F("hTTbar","",gIndex, -0.5, gIndex-0.5);
  TH1F* hTTW   = new TH1F("hTTW"  ,"",gIndex, -0.5, gIndex-0.5);
  TH1F* hTTH   = new TH1F("hTTH"  ,"",gIndex, -0.5, gIndex-0.5);
  THStack* mc  = new THStack("mc","mc");
  vector<Point>::iterator point;
  cout << "TTbar size is " << fTTbarMC.size() << endl;
  for (point = fTTbarMC.begin(); point != fTTbarMC.end(); ++point)
    hTTbar->Fill( mainCluster.FindUnclusterizableCluster( *point), 37000*point->fW);
  for (point = fTTWMC.begin(); point != fTTWMC.end(); ++point)
    hTTW->Fill( mainCluster.FindUnclusterizableCluster( *point), 37000*point->fW);
  for (point = fTTHMC.begin(); point != fTTHMC.end(); ++point)
    hTTH->Fill( mainCluster.FindUnclusterizableCluster( *point), 37000*point->fW);
  cout << hTTbar->Integral() << " " << hTTbar->GetEntries() << endl;
  hTTbar->SetFillColor( kRed     );
  hTTH->SetFillColor( kBlue    );
  hTTW->SetFillColor( kMagenta );

  mc->Add( hTTbar ); mc->Add(hTTW); mc->Add(hTTH);
  mc->Draw("HIST");

  for (int k = 0; k < gIndex; ++k){
    cout << hTTbar->GetBinContent(k+1) << "  " 
	 << hTTW  ->GetBinContent(k+1) << "  " 
	 << hTTH  ->GetBinContent(k+1) << endl;
  }

#ifdef SIGNIFICANCE_H
  Significance c;
  c.Test();
#endif 
  cout << "it works" << endl;
  
  vector<TH1*> bkgs;
  bkgs.push_back(hTTbar);
  bkgs.push_back(hTTW  );
  //  MakeSimpleCard card(hTTH, bkgs, "datacard_newBinning", 37000., false);
  MakeSimpleCard card(hTTH, bkgs, "datacard_recursiveclustering", 1., false);
  card.doCard();

  return;
}

void RecursiveClustering::StoreToFile()
{
  TFile* binning = TFile::Open("binning.root","recreate");
  TH2F*  hBinning = new TH2F("hBinning","",1000,-1.,1.,1000,-1.,1.);
  for (Int_t binx = 0; binx < hBinning->GetXaxis()->GetNbins(); ++binx){
      for (Int_t biny = 0; biny < hBinning->GetYaxis()->GetNbins(); ++biny){
	Double_t x  = hBinning->GetXaxis()->GetBinCenter(binx);
	Double_t y  = hBinning->GetYaxis()->GetBinCenter(biny);
	Int_t bin = hBinning->GetBin(binx,biny);
	hBinning->SetBinContent(bin, mainCluster.FindUnclusterizableCluster(Point(x,y,-1)));
      }
  }
  hBinning->Write();
  binning->Close();


}


void RecursiveClustering::VoronoiPlot()
{
  vector<TGraph*> graphs; graphs.clear();
  vector<Double_t>* X = new vector<Double_t>[gIndex];
  vector<Double_t>* Y = new vector<Double_t>[gIndex]; 

  cout << "Calculating points" << endl;
  for (Double_t x = -1; x < 1.; x = x + 1e-3){
    
      for (Double_t y = -1; y < 1.; y = y + 1e-3){
	Int_t k = mainCluster.FindUnclusterizableCluster(Point(x,y,-1));
	X[k].push_back(x);
	Y[k].push_back(y);	
      }
  }



  TCanvas* c = new TCanvas();
  c->cd();
  gStyle->SetOptStat(0);
  TH1F* hDummy = new TH1F("hDummy","",2,-1,1);
  hDummy->SetBinContent(1, 1.);
  hDummy->SetBinContent(2,-1.);
  hDummy->SetLineColor(kWhite);
  hDummy->GetYaxis()->SetRangeUser(-1.,1.);
  hDummy->Draw();
  cout << "Done... now plotting" << endl;
  cout << fCentroids.size() << endl;

  for (Int_t k = 0; k < gIndex; ++k){
    graphs.push_back(new TGraph( X[k].size(), &X[k][0], &Y[k][0] ));
    graphs[k]->SetMarkerColor(k+1);    
    graphs[k]->SetMarkerStyle(6);
    graphs[k]->Draw("PSAME");
  }  

}



void RecursiveClustering::VoronoiPlot2(Int_t level)
{
  vector<TGraph*> graphs; graphs.clear();
  vector<Double_t>* X = new vector<Double_t>[gIndex];
  vector<Double_t>* Y = new vector<Double_t>[gIndex]; 

  vector<TString> stringList; stringList.clear();
  for (Double_t x = -1; x < 1.; x = x + 1e-3){
      for (Double_t y = -1; y < 1.; y = y + 1e-3){
	TString nam = mainCluster.FindSubClusterName(Point(x,y,-1), level);
	vector<TString>::iterator itr = std::find( stringList.begin(), stringList.end(), nam);
	if ( itr == stringList.end()){
	  stringList.push_back(nam);
	}
	itr = std::find( stringList.begin(), stringList.end(), nam);
	X[itr-stringList.begin()].push_back(x);
	Y[itr-stringList.begin()].push_back(y);	
      }
  }


  TCanvas* c = new TCanvas();
  c->cd();
  gStyle->SetOptStat(0);
  TH1F* hDummy = new TH1F("hDummy","",2,-1,1);
  hDummy->SetBinContent(1, 1.);
  hDummy->SetBinContent(2,-1.);
  hDummy->SetLineColor(kWhite);
  hDummy->GetYaxis()->SetRangeUser(-1.,1.);
  hDummy->Draw();


  for (Int_t k = 0; k < gIndex; ++k){
    graphs.push_back(new TGraph( X[k].size(), &X[k][0], &Y[k][0] ));
    graphs[k]->SetMarkerColor(k+1);    
    graphs[k]->SetMarkerStyle(6);
    graphs[k]->Draw("PSAME");
  }  

}
