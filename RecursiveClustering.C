#include "RecursiveClustering.h"
#include "MakeSimpleCard.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include <sstream>
#include <TString.h> 

using namespace std;

bool SortX( Point i, Point j)
{
  return ( i.fX < j.fX );
}
bool SortY( Point i, Point j)
{
  return ( i.fY < j.fY );
}

Int_t gIndex = 0;

Cluster::Cluster( vector<Point> TTbar, vector<Point> TTH  , vector<Point> TTW, Int_t k )
{
  cout << "Initialising cluster with " << TTH.size() 
       << " " << TTW.size() << " " << TTbar.size() << endl;
  fTTbar  = TTbar;
  fTTH    = TTH  ;	
  fTTW    = TTW  ;	

  fData = fTTbar; 
  fData.insert( fData.end(), fTTW.begin(),fTTW.end()); 
  fData.insert( fData.end(), fTTH.begin(),fTTH.end());

  fIsClusterizable = true;
  fK = k;


}
void Cluster::MakeAssignment()
{
  fTgt.clear();
  // Assigns each data point to its cluster :)
  for (int n = 0; n < fData.size(); ++n){
    fTgt.push_back(FindSubCluster( fData[n] ));
  }
}

Int_t Cluster::FindSubCluster( Point point )
{
  Int_t cluster = -1;
  Double_t dist = 10000;
  for (int k = 0; k < fCentroids.size(); ++k){
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
  if (fIsClusterizable){
    return SubClusters[FindSubCluster(point)].FindUnclusterizableCluster(point);
  }
  else{
    return fIndex;
  }

}

Int_t Cluster::CalculateCentroids()
{
  Int_t theyChange = 1;
  vector<Double_t> fMWght;
  vector<Point> fOldCentroids = fCentroids;

  for (int k = 0; k < fCentroids.size(); ++k){
    fCentroids[k].fX = 0;
    fCentroids[k].fY = 0;
    fCentroids[k].fW = -1000000; // weights should never be used for centroids
    fMWght.push_back(0);
  }

  for (int n = 0; n < fData.size(); ++n){
    fCentroids[fTgt[n]].fX    +=  fData[n].fW * fData[n].fX;
    fCentroids[fTgt[n]].fY    +=  fData[n].fW * fData[n].fY;
    fMWght[fTgt[n]] +=  fData[n].fW;
  }

  Double_t di = 0;
  
  for (int k = 0; k < fCentroids.size(); ++k){
    fCentroids[k].fX /= fMWght[k];
    fCentroids[k].fY /= fMWght[k];
    di += d(fCentroids[k].fX,fCentroids[k].fY, fOldCentroids[k].fX, fOldCentroids[k].fY);
  }
  if (di < 1.0e-8) theyChange =  -1;
  // cout << "[I]: New centroids are (" << fMx[0] << "," << fMy[0]
  //      << ") and (" << fMx[1] << "," << fMy[1] << ")" << endl;


  return theyChange;

}


void Cluster::recluster()
{
  TRandom3* r = new TRandom3();
  Double_t maxX = (*max_element(fData.begin(), fData.end(),SortX)).fX;
  Double_t minX = (*min_element(fData.begin(), fData.end(),SortX)).fX;
  Double_t maxY = (*max_element(fData.begin(), fData.end(),SortY)).fY;
  Double_t minY = (*min_element(fData.begin(), fData.end(),SortY)).fY;

  // Random inizialitation of centroids
  for (int k = 0; k < fK; ++k){
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
  
  
  cout << "Classificating points and building subclusters" << endl;
  for (int k = 0; k < fK; ++k){
    cout << "Cluster " << k << endl;
    vector<Point> TTbar;
    vector<Point> TTW;
    vector<Point> TTH;
    for (int n = 0; n < fTTH.size(); ++n){
      if ( FindSubCluster( fTTH[n] ) == k) TTH.push_back( fTTH[n]);
    }
    for (int n = 0; n < fTTW.size(); ++n){
      if ( FindSubCluster( fTTW[n] ) == k) TTW.push_back( fTTH[n]);
    }
    for (int n = 0; n < fTTbar.size(); ++n){
      if ( FindSubCluster( fTTbar[n] ) == k) TTbar.push_back( fTTH[n]);
    }
    cout << TTH.size() << " " << TTbar.size() << " " << TTW.size() << endl;
    if (TTH.size() == 0) fIsClusterizable = false;
    else if (TTbar.size() == 0) fIsClusterizable = false;
    else if (TTW  .size() == 0) fIsClusterizable = false;
    else{
      if (36400*TTH.size()*TTH[0].fW < 5.)
	fIsClusterizable = false;
    }
    if (!fIsClusterizable){
      cout << "Not clusterizable anymore " << endl;
      cout << "This will be subcluster " << gIndex <<  endl;
      SubClusters.clear();
      fIndex = gIndex;
      gIndex++;
      break;
    }

      
    Cluster subCluster = Cluster( TTbar,  TTH  ,  TTW, fK);
    subCluster.recluster();
    SubClusters.push_back(subCluster);
  }


  delete r;

}


RecursiveClustering::RecursiveClustering(Int_t k)
{
  fK = k;
  readFromFiles();
  StartTheThing();
  
  cout << "Produced " << gIndex << " clusters" << endl;
}

void RecursiveClustering::StartTheThing()
{
  mainCluster = Cluster(fTTbar, fTTH, fTTW, fK);
  mainCluster.recluster();
}

void RecursiveClustering::readFromFiles()
{
  
  fTTbar.clear(); fTTH.clear(); fTTW.clear();
  ifstream f; f.open("data/ttbar2.txt");
  Int_t count = 0;
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if ( TMath::Abs(x) > 80. ) continue;
    if ( TMath::Abs(y) > 80. ) continue;
    Point point = Point(x,y,w);
    if (count%2 == 0)
      fTTbar.push_back(point);
    else
      fTTbarMC.push_back(point);
    count++;
  }
  f.close();
  
  f.open("data/tth.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    Point point = Point(x,y,w);
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
    Point point = Point(x,y,w);
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

  for (point = fTTbarMC.begin(); point != fTTbarMC.end(); ++point)
    hTTbar->Fill( mainCluster.FindUnclusterizableCluster( *point), point->fW);
  for (point = fTTWMC.begin(); point != fTTWMC.end(); ++point)
    hTTW->Fill( mainCluster.FindUnclusterizableCluster( *point), point->fW);
  for (point = fTTHMC.begin(); point != fTTHMC.end(); ++point)
    hTTH->Fill( mainCluster.FindUnclusterizableCluster( *point), point->fW);
  hTTbar->SetFillColor( kRed     );
  hTTH->SetFillColor( kBlue    );
  hTTW->SetFillColor( kMagenta );

  mc->Add( hTTbar ); mc->Add(hTTW); mc->Add(hTTH);
  mc->Draw("HIST");

  vector<TH1*> bkgs;
  bkgs.push_back(hTTbar);
  bkgs.push_back(hTTW  );

  MakeSimpleCard card(hTTH, bkgs, "datacard", true);
  card.doCard();

  return;
}


void RecursiveClustering::VoronoiPlot()
{
  vector<TGraph*> graphs; graphs.clear();
  vector<Double_t>* X = new vector<Double_t>[gIndex];
  vector<Double_t>* Y = new vector<Double_t>[gIndex]; 


  for (Double_t x = -1; x < 1.; x = x + 1e-3){
      for (Double_t y = -1; y < 1.; y = y + 1e-3){
	Int_t k = mainCluster.FindUnclusterizableCluster(Point(x,y,-1));
	X[k].push_back(x);
	Y[k].push_back(y);	
      }
  }


  TCanvas* c = new TCanvas();
  c->cd();
  TH1F* hDummy = new TH1F("hDummy","",2,-1,1);
  hDummy->SetBinContent(1, 1.);
  hDummy->SetBinContent(2,-1.);
  hDummy->SetLineColor(kWhite);
  hDummy->Draw();

  for (Int_t k = 0; k < gIndex; ++k){
    graphs.push_back(new TGraph( X[k].size(), &X[k][0], &Y[k][0] ));
    graphs[k]->SetMarkerColor(k+1);    
    graphs[k]->SetMarkerStyle(6);
    graphs[k]->Draw("PSAME");
  }  

}
