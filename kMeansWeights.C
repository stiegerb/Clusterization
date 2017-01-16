#include "kMeansWeights.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>

using namespace std;

kMeansWeights::kMeansWeights(Int_t k)
{
  fK = k;
  ostringstream convert;   // stream used for the conversion
  convert << fK;      // insert the textual representation of 'Number' in the characters in the stream
  fKname = TString(convert.str()); // set 'Result' to the contents of the stream

  readFromFiles();

}

void kMeansWeights::readFromFiles()
{
  
  fN = 0;
  ifstream f; f.open("data/ttbartest.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    fX.push_back(x);
    fY.push_back(y);
    fW.push_back(w);
    fN ++;
  }
  f.close();


}

// Double_t kMeansWeights::getMax( Double_t* data, Double_t coeff)
// {
//   Double_t max = -9999;
//   for (int n = 0; n < fN; ++n){
//     if ( coeff*data[n] > max)
//       max = data[n];
//   }
//   if ( max == -9999) cout << "[E]: Maximum hasnt been calculated" << endl;
//   return max;
// }

Double_t kMeansWeights::d(Double_t x, Double_t y, Double_t m_x, Double_t m_y) 
{
  return TMath::Sqrt( (x-m_x) * (x-m_x) + (y-m_y) * (y-m_y) );
}


void kMeansWeights::MakeAssignment()
{
  // Assigns each data point to its cluster :)
  for (int n = 0; n < fN; ++n){
    Int_t cluster = -1;
    Double_t dist = 10000;
    for (int k = 0; k < fK; ++k){
      Double_t dst = d( fX[n], fY[n], fMx[k], fMy[k]);
      //      cout << "A" << fX[k] << " " <<  fY[k] << " " <<  fMx[k] << " " <<  fMy[k] << endl;
      // cout << "B" << k << " " << dist<< " " << dst << endl;
      if (dst < dist){
	cluster = k;
	dist = dst;
      }
    }
    if (cluster < 0){ cout << "[E]: Nearest cluster not found" << endl;}
    fTgt[n] = cluster; 
  }
}

Int_t kMeansWeights::CalculateCentroids()
{
  Int_t theyChange = 1;
  Double_t* fMWght = new Double_t[fK];
  Double_t* fOldMx = new Double_t[fK];
  Double_t* fOldMy = new Double_t[fK];

  for (int k = 0; k < fK; ++k){
    fOldMx[k] = fMx[k];
    fOldMy[k] = fMy[k];
    fMx[k] = 0; fMy[k] = 0;
    fMWght[k] = 0;
  }

  for (int n = 0; n < fN; ++n){
    fMx[fTgt[n]]    +=  fW[n] * fX[n];
    fMy[fTgt[n]]    +=  fW[n] * fY[n];
    fMWght[fTgt[n]] +=  fW[n];
  }

  Double_t di = 0;
  
  for (int k = 0; k < fK; ++k){
    fMx[k] /= fMWght[k];
    fMy[k] /= fMWght[k];
    di += d(fMx[k],fMy[k], fOldMx[k], fOldMy[k]);
  }
  if (di < 1.0e-8) theyChange =  -1;
  cout << "[I]: New centroids are (" << fMx[0] << "," << fMy[0]
       << ") and (" << fMx[1] << "," << fMy[1] << ")" << endl;

  delete fMWght;
  delete fOldMx;
  delete fOldMy;

  return theyChange;

}

void kMeansWeights::Init()
{
  fMx  = new Double_t[fK];
  fMy  = new Double_t[fK];
  // fX   = new Double_t[fN];
  // fY   = new Double_t[fN];
  // fW   = new Double_t[fN];
  fTgt = new Int_t[fN];

}

void kMeansWeights::fit()
{
  TRandom3* r = new TRandom3();
  Double_t maxX = (*max_element(fX.begin(), fX.end()));
  Double_t minX = (*min_element(fX.begin(), fX.end()));
  Double_t maxY = (*max_element(fY.begin(), fY.end()));
  Double_t minY = (*min_element(fY.begin(), fY.end()));

  Init();

  // Random inizialitation of centroids
  for (int k = 0; k < fK; ++k){
    fMx[k] = r->Uniform(minX, maxX);
    fMy[k] = r->Uniform(minY, maxY);
  }
  // Cluster assignation
  Int_t maxIt = 999999;
  Int_t it = 0;

  while (it < maxIt){
    it++;
    MakeAssignment();
    if (CalculateCentroids() < 0) break;

  }
  if ( maxIt == it) cout << "[W]: Maximum number of iterations reached!! ";
  
  fDataX = new vector<Double_t>[fK];
  fDataY = new vector<Double_t>[fK];
  for (int n = 0; n < fN; ++n){
    fDataX[fTgt[n]].push_back( fX[n] );
    fDataY[fTgt[n]].push_back( fY[n] );
  }

  delete r;

}

void kMeansWeights::ScatterPlot()
{

  TCanvas* c = new TCanvas("kmeans"+fKname, "kmeans"+fKname, 800, 800);
  c->cd();
  std::vector<TGraph*> graphs; graphs.clear();
  for(int igr=0; igr<fK; ++igr)
    {
      graphs.push_back(new TGraph( fDataX[igr].size(), &fDataX[igr][0], &fDataY[igr][0] ));
      graphs[igr]->SetMarkerColor(igr+1);    
      graphs[igr]->SetMarkerStyle(kFullCircle);
      
      igr==0 ? graphs[igr]->Draw("AP") : graphs[igr]->Draw("PSAME");
      
    }
  c->Print("kmeans"+fKname+".png");
  c->Print("kmeans"+fKname+".pdf");

}
