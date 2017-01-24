#include "kMeansWeights.h"
//#ifdef MAKESIMPLECARD_H
#include "MakeSimpleCard.h"
//#endif
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
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <sstream>

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
  ifstream f; f.open("data/TTSingleLeptonMarco.txt");
  Int_t count = 0;
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if (count%2 == 0){
      fX.push_back(x);
      fY.push_back(y);
      fW.push_back(2*w);
      cout << x << " "  << y << " " << w  << endl;
      fN ++;
    }
    else{
      fTTbarX.push_back(x);
      fTTbarY.push_back(y);
      fTTbarW.push_back(2*w);
    }
    count++;
  }
  f.close();
  
  f.open("data/ttH.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if (count%2 == 0){
      fX.push_back(x);
      fY.push_back(y);
      fW.push_back(2*w);
      fN ++;
    }
    else{
      fTTHX.push_back(x);
      fTTHY.push_back(y);
      fTTHW.push_back(2*w);
    }
    count++;
  }
  f.close();
  
  f.open("data/ttw.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if (count%2 == 0){
      fX.push_back(x);
      fY.push_back(y);
      fW.push_back(2*w);
      fN ++;
    }
    else{
      fTTWX.push_back(x);
      fTTWY.push_back(y);
      fTTWW.push_back(2*w);
    }
    count++;
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
  //return TMath::Abs(x-m_x) + TMath::Abs(y-m_y);
}


void kMeansWeights::MakeAssignment()
{
  // Assigns each data point to its cluster :)
  for (int n = 0; n < fN; ++n){
    Int_t cluster = -1;
    Double_t dist = 10000;
    for (int k = 0; k < fK; ++k){
      Double_t dst = d( fX[n], fY[n], fMx[k], fMy[k]);
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
    cout << fMx[k] << " " << fMy[k] << endl; 
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
  TH1F* hDummy = new TH1F("hDummy","",2,-1,1);
  hDummy->SetBinContent(1, 1.);
  hDummy->SetBinContent(2,-1.);
  hDummy->SetLineColor(kWhite);
  hDummy->Draw();
  
  std::vector<TGraph*> graphs; graphs.clear();
  for(int igr=0; igr<fK; ++igr){
      graphs.push_back(new TGraph( fDataX[igr].size(), &fDataX[igr][0], &fDataY[igr][0] ));
      graphs[igr]->SetMarkerColor(igr+1);    
      graphs[igr]->SetMarkerStyle(kFullCircle);
      
      //igr==0 ? graphs[igr]->Draw("AP") : graphs[igr]->Draw("PSAME");
      graphs[igr]->Draw("PSAME");
      
    }
  graphs[0]->GetXaxis()->SetRangeUser(-1.,1.);
  graphs[0]->GetYaxis()->SetRangeUser(-1.,1.);
  c->Print("kmeans"+fKname+".png");
  c->Print("kmeans"+fKname+".pdf");

}

Int_t kMeansWeights::GetCluster(Double_t x, Double_t y)
{
  Int_t cluster = -1;
  Double_t dist = 99999;
  for (int k = 0; k < fK; ++k){
    Double_t dst = d(x, y, fMx[k], fMy[k]);
    if (dst < dist){
      cluster = k;
      dist = dst;
    }
  }
  if (cluster < 0){ cout << "[E]: Nearest cluster not found" << endl;}
  return cluster;
}
Int_t kMeansWeights::classicalBinning(Double_t x, Double_t y){
  if ((-1. < x) && (x <= 1) && (-1.0 < y) && (y <= -0.2))  return 0;
  else if ((-1. < x ) && ( x <= 1)   && (-0.2 < y) && ( y  <=  0.1))  return 1;
  else if ((-1. < x ) && ( x <= 0.3) && (0.1  < y) && ( y  <=  0.4))  return 2;
  else if ((0.3 < x ) && ( x <= 1.)  && (0.1  < y) && ( y  <=  0.4))  return 3;
  else if ((-1. < x ) && ( x <= 0.1) && (0.4  < y) && ( y  <=  1.0))  return 4;
  else if ((0.1 < x ) && ( x <= 0.4) && (0.4  < y) && ( y  <=  1.0))  return 5;
  else if ((0.4 < x ) && ( x <= 1.)  && (0.4  < y) && ( y  <=  1.0))  return 6;
  else {
    cout <<  "one bin is missing " << x << " " << y << endl;
    return -1;
  }
}



void kMeansWeights::makeHistos()
{
  TH1F* hTTb = new TH1F("hTTb","", fK, -0.5, fK-0.5);
  TH1F* hTTH = new TH1F("hTTH","", fK, -0.5, fK-0.5);
  TH1F* hTTW = new TH1F("hTTW","", fK, -0.5, fK-0.5);

  TH1F* hTTbOld = new TH1F("hTTbOld","", fK, -0.5, fK-0.5);
  TH1F* hTTHOld = new TH1F("hTTHOld","", fK, -0.5, fK-0.5);
  TH1F* hTTWOld = new TH1F("hTTWOld","", fK, -0.5, fK-0.5);

  for (size_t n = 0; n < fTTbarX.size(); ++n){
    hTTb->Fill( GetCluster(fTTbarX[n], fTTbarY[n]), fTTbarW[n]);
    hTTbOld->Fill( classicalBinning(fTTbarX[n], fTTbarY[n]), fTTbarW[n]);
  }
  cout << hTTbOld->Integral() << " " << hTTbOld->GetEntries() << endl;
  for (size_t n = 0; n < fTTHX.size(); ++n){
    hTTH->Fill( GetCluster(fTTHX[n], fTTHY[n]), fTTHW[n]);
    hTTHOld->Fill( classicalBinning(fTTHX[n], fTTHY[n]), fTTHW[n]);
    
  }
  for (size_t n = 0; n < fTTWX.size(); ++n){
    hTTW->Fill( GetCluster(fTTWX[n], fTTWY[n]), fTTWW[n]);
    hTTWOld->Fill( classicalBinning(fTTWX[n], fTTWY[n]), fTTWW[n]);
  }
  THStack* mc = new THStack();
  THStack* mcOld = new THStack();
  hTTb->SetFillColor( kRed     );
  hTTH->SetFillColor( kBlue    );
  hTTW->SetFillColor( kMagenta );
  hTTbOld->SetFillColor( kRed     );
  hTTHOld->SetFillColor( kBlue    );
  hTTWOld->SetFillColor( kMagenta );
  mc->Add( hTTH);  mc->Add( hTTb);  mc->Add( hTTW);
  mcOld->Add( hTTHOld);  mcOld->Add( hTTbOld);  mcOld->Add( hTTWOld);
  TCanvas* c1 = new TCanvas();
  mc->Draw("HIST");
  TCanvas* c2 = new TCanvas();
  mcOld->Draw("HIST");

  TH1F* hTTb2 = (TH1F*) hTTb->Clone("hTTb2");
  TH1F* hTTW2 = (TH1F*) hTTW->Clone("hTTW2");
  TH1F* hTTH2 = (TH1F*) hTTH->Clone("hTTH2");
  TH1F* hTTbOld2 = (TH1F*) hTTbOld->Clone("hTTbOld2");
  TH1F*	hTTHOld2 = (TH1F*) hTTHOld->Clone("hTTHOld2");
  TH1F*	hTTWOld2 = (TH1F*) hTTWOld->Clone("hTTWOld2");
  
  // hTTb2->Add(hTTW2); 
  // cout << "New " << hTTH2->KolmogorovTest(hTTb2,"M") << endl;
  // hTTbOld2->Add(hTTWOld2);
  // cout << "Old " << hTTHOld2->KolmogorovTest(hTTbOld2,"M") << endl;
    
  //#ifdef MAKESIMPLECARD_H
  vector<TH1*> bkgs;
  bkgs.push_back(hTTbOld2);
  bkgs.push_back(hTTWOld2);
  MakeSimpleCard card(hTTHOld2, bkgs, "datacard_classical", 37000., false);
  card.doCard();

  vector<TH1*> bkgs_c;
  bkgs_c.push_back(hTTb2);
  bkgs_c.push_back(hTTW2);
  MakeSimpleCard card_c(hTTH2, bkgs_c, "datacard_kmeansweights", 37000., false);
  card_c.doCard();
  //#endif

  #ifdef SIGNIFICANCE_H
  Significance c;
  c.Test();
#endif 

  return;
}


void kMeansWeights::MCScatterPlots()
{
  TGraph* gr1 = new TGraph( fTTbarX.size(), &fTTbarX[0], &fTTbarY[0]);
  TGraph* gr2 = new TGraph( fTTHX.size()  , &fTTHX[0]  , &fTTHY[0]  );
  TGraph* gr3 = new TGraph( fTTWX.size()  , &fTTWX[0]  , &fTTWY[0]  );

  gr1->SetMarkerStyle( 6 );
  gr2->SetMarkerStyle( 6 );
  gr3->SetMarkerStyle( 6 );


  TCanvas* c1 = new TCanvas(); gr1->Draw("A,P");
  TCanvas* c2 = new TCanvas(); gr2->Draw("A,P");
  TCanvas* c3 = new TCanvas(); gr3->Draw("A,P");
  gr1->GetXaxis()->SetRangeUser(-1.,1.);
  gr1->GetYaxis()->SetRangeUser(-1.,1.);
  gr2->GetXaxis()->SetRangeUser(-1.,1.);
  gr2->GetYaxis()->SetRangeUser(-1.,1.);
  gr3->GetXaxis()->SetRangeUser(-1.,1.);
  gr3->GetYaxis()->SetRangeUser(-1.,1.);

  
  
  return;
}


void kMeansWeights::VoronoiPlot()
{
  vector<TGraph*> graphs; graphs.clear();
  vector<Double_t>* X = new vector<Double_t>[fK];
  vector<Double_t>* Y = new vector<Double_t>[fK]; 

  cout << "Calculating points" << endl;
  for (Double_t x = -1; x < 1.; x = x + 1e-3){
    
      for (Double_t y = -1; y < 1.; y = y + 1e-3){
	Int_t k = GetCluster(x,y);
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

  for (Int_t k = 0; k < fK; ++k){
    graphs.push_back(new TGraph( X[k].size(), &X[k][0], &Y[k][0] ));
    graphs[k]->SetMarkerColor(k+1);    
    graphs[k]->SetMarkerStyle(6);
    graphs[k]->Draw("PSAME");
  }  

}
