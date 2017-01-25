#include "TargettedClustering.h"
#include "MakeSimpleCard.h"
#include <TFile.h>
//#ifdef SIGNIFICANCE_H
#include "Significance.h"
//#endif
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
#include "TSystem.h"
#include "TTree.h"
#include "TFormula.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TROOT.h"

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

Cluster::Cluster( vector<Point> TTbar, vector<Point> TTH  , vector<Point> TTW, Int_t k , TString name)
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
  return ( fFormula.Eval( point.fX, point.fY) > 0.) ? 0 : 1;
  
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

double Cluster::FOMforClusterCut( const double *par)
{
  fFormula.SetParameters(par);
  Double_t s1 = 0.;
  Double_t s2 = 0.;
  Double_t ttbar1 = 0.;
  Double_t ttbar2 = 0.;
  Double_t ttw1 = 0.;
  Double_t ttw2 = 0.;
  vector<Point>::iterator point;
  for ( point = fTTbar.begin(); point != fTTbar.end(); ++point){
    if ( fFormula.Eval( point->fX, point->fY ) > 0)
      ttbar1 += point->fW;
    else
      ttbar2 += point->fW;
  }
  for ( point = fTTW.begin(); point != fTTW.end(); ++point){
    if ( fFormula.Eval( point->fX, point->fY ) > 0)
      ttw1 += point->fW;
    else
      ttw2 += point->fW;
  }
  for ( point = fTTH.begin(); point != fTTH.end(); ++point){
    if ( fFormula.Eval( point->fX, point->fY ) > 0)
      s1 += point->fW;
    else
      s2 += point->fW;
  }

  ttbar1 *= 37000.;
  ttbar2 *= 37000.;
  ttw1   *= 37000.;
  ttw2   *= 37000.;
  s1 *= 37000.;
  s2 *= 37000.;

  Double_t pvalue1 = TMath::Poisson( s1+ttbar1+ttw1, ttbar1+ttw1 );
  Double_t pvalue2 = TMath::Poisson( s2+ttbar2+ttw2, ttbar2+ttw2 );
  
  double result = 0.;
  if ( (s1 * s2 * ttbar1 * ttbar2 * ttw1 * ttw2) == 0) result = 99999;
  else result = pvalue1 * pvalue2;
  return result;
}

void Cluster::recluster()
{
  cout << "Reclustering " << endl;
  TRandom3* r = new TRandom3();


  fFormula = TFormula(fName+"formula", "[0] + [1]*x < y");
  fFormula.SetParameter(0,0.);
  fFormula.SetParameter(1,1.);

  vector<Point> subCentroidList; subCentroidList.clear();
  vector<double> significancesList; significancesList.clear();
  vector<double> partialSignificancesList; partialSignificancesList.clear();
  
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit","Combined");
  ROOT::Math::Functor functor(this, &Cluster::FOMforClusterCut, 2);
  min->SetFunction(functor);

  min->SetVariable(0,"a",r->Uniform(-1.,1.),0.01);
  min->SetVariable(1,"b",r->Uniform(-5, 5),0.01);
  min->Minimize();
  const double *xs = min->X();
  
  while ( FOMforClusterCut(xs) > 2.){
    min->SetVariable(0,"a",r->Uniform(-1.,1.),0.01);
    min->SetVariable(1,"b",r->Uniform(-20.,20.),0.01);
    min->Minimize();
    const double *xs = min->X();
  }
  
  cout << "8==================D" << endl;
  cout << "Optimal cuts are..." << endl;
  cout << xs[0] << " " << xs[1] << endl;
  cout << "8==================D" << endl;
  

  

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

    vector<double> yields;
    yields.push_back(36400*TTH.size()*TTH[0].fW);
    yields.push_back(36400*(TTW.size()*TTW[0].fW + TTbar.size()*TTbar[0].fW));
    Significance c(yields);
    cout << "The significance is " << c.getSignificance("pvalue") << endl;
    if (!fIsClusterizable){
      cout << "Not clusterizable anymore " << endl;
      cout << "This will be subcluster " << gIndex <<  endl;
      cout << "The name is " << fName << endl;
      SubClusters.clear();
      fIndex = gIndex;
      gIndex++;
      significancesList.push_back(c.getSignificance("pvalue"));
      break;
    }
    else{
      cout << "Apparently subcluster " << k << " from " << fName 
	   << " is huge!!! (thats what she said), so theres not showstopper not to keep clustering" << endl;
      Cluster subCluster = Cluster( TTbar,  TTH  ,  TTW, fK, fName + Form("%u",k));
      cout << "Subcluster is produced " << endl;
      SubClusters.push_back(subCluster);
    }
  }

  for (size_t k = 0; k < SubClusters.size(); ++k){
      cout << "Cluster " << fName << " is huge!!! (thats what she said)" << endl;
      SubClusters[k].recluster();
  }


  delete r;
  return;

}


TargettedClustering::TargettedClustering(Int_t k, Int_t nLep):
  fK(k),
  nLep_(nLep)
{
  readFromFiles();
  StartTheThing();
  
  cout << "Produced " << gIndex << " clusters" << endl;
  Point point(0.353811, 0.456623,-1.);
  cout << "Point is " << mainCluster.FindUnclusterizableCluster(point) <<endl;
}

void TargettedClustering::StartTheThing()
{
  gROOT->LoadMacro("Significance.C+");
  gROOT->LoadMacro("MakeSimpleCard.C+");

  mainCluster = Cluster(fTTbar, fTTH, fTTW, fK, "A");
  mainCluster.recluster();
  cout << "Final list of significances: " << endl;
  double combinedSignificance(1.);
  for(auto& significance : fSignificances)
    {
      combinedSignificance*=significance;
      cout << significance << endl;
    }
  cout << "Final combined significance: " << combinedSignificance << endl;
  StoreToFile();
}

void TargettedClustering::readFromFiles()
{
  
  fTTbar.clear(); fTTH.clear(); fTTW.clear();
  ifstream f;
  nLep_==3 ? f.open("data/ttbar3l.txt") : f.open("data/TTSingleLeptonMarco.txt");
  Int_t count = 0;
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if ( TMath::Abs(x) > 2.) continue;
    if ( TMath::Abs(y) > 2.) continue;
    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTbar.push_back(point);
    else
      fTTbarMC.push_back(point);
    count++;
  }
  f.close();
  cout << "TTbar events " << fTTbar.size() << endl;
  nLep_==3 ? f.open("data/tth3l.txt") : f.open("data/ttH.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if (w == 0) continue;
    if ( TMath::Abs(x) > 2.) continue;
    if ( TMath::Abs(y) > 2.) continue;

    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTH.push_back(point);
    else
      fTTHMC.push_back(point);
    count++;
  }
  f.close();
  cout << "TTH events " << fTTH.size() << endl;
  nLep_==3 ? f.open("data/ttw3l.txt") : f.open("data/ttw.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> x >> y >> w;
    if (w == 0) continue;
    if ( TMath::Abs(x) > 2.) continue;
    if ( TMath::Abs(y) > 2.) continue;

    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTW.push_back(point);
    else
      fTTWMC.push_back(point);
    count++;
  }
  f.close();

}

// // Double_t TargettedClustering::getMax( Double_t* data, Double_t coeff)
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

void TargettedClustering::Test()
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

//#ifdef SIGNIFICANCE_H
//  cout << "Now significance test" << endl;
//  Significance c;
//  c.Test();
//#endif 
//
//#ifdef MAKESIMPLECARD_H  
//  cout << "Now simple card test " << endl;
//  vector<TH1*> bkgs;
//  bkgs.push_back(hTTbar);
//  bkgs.push_back(hTTW  );
//
//  MakeSimpleCard card(hTTH, bkgs, "datacard_recursiveclustering", 1., false);
//  card.doCard();
//#endif

}

void TargettedClustering::StoreToFile()
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


  // Final significance per bin
  double indexes[fSignificances.size()];
  double signifs[fSignificances.size()];
  for(int idx=0; idx<fSignificances.size(); ++idx)
    {
      indexes[idx]=idx;
      signifs[idx]=fSignificances[idx];
    }
  TGraph* g = new TGraph(fSignificances.size(), indexes, signifs);
  g->SetName("significancePerBin_finalSet");
  g->Write();

  binning->Close();
}


void TargettedClustering::VoronoiPlot()
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


  for (Int_t k = 0; k < gIndex; ++k){
    graphs.push_back(new TGraph( X[k].size(), &X[k][0], &Y[k][0] ));
    graphs[k]->SetMarkerColor(k+1);    
    graphs[k]->SetMarkerStyle(6);
    graphs[k]->Draw("PSAME");
  }  

}



void TargettedClustering::VoronoiPlot2(Int_t level)
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



Double_t TargettedClustering::SignificanceAtLevel(Int_t level)
{
  vector<TGraph*> graphs; graphs.clear();
  vector<TString> stringList; stringList.clear();
  vector<Point>::iterator point; 
  TH1F* hTTbar = new TH1F(Form("hTTbarLevel%d",level) , "", gIndex, -0.5, gIndex-0.5);
  TH1F* hTTH   = new TH1F(Form("hTTHLevel%d",level)   , "", gIndex, -0.5, gIndex-0.5);
  TH1F* hTTW   = new TH1F(Form("hTTWLevel%d",level)   , "", gIndex, -0.5, gIndex-0.5);

  for (point = fTTbar.begin(); point != fTTbar.end(); ++point){
	TString nam = mainCluster.FindSubClusterName(*point, level);
	vector<TString>::iterator itr = std::find( stringList.begin(), stringList.end(), nam);
	if ( itr == stringList.end()){
	  stringList.push_back(nam);
	}
	itr = std::find( stringList.begin(), stringList.end(), nam);
	hTTbar->Fill(itr-stringList.begin(), point->fW);
  }

  for (point = fTTW.begin(); point != fTTW.end(); ++point){
	TString nam = mainCluster.FindSubClusterName(*point, level);
	vector<TString>::iterator itr = std::find( stringList.begin(), stringList.end(), nam);
	if ( itr == stringList.end()){
	  stringList.push_back(nam);
	}
	itr = std::find( stringList.begin(), stringList.end(), nam);
	hTTW->Fill(itr-stringList.begin(), point->fW);
  }

  for (point = fTTH.begin(); point != fTTH.end(); ++point){
	TString nam = mainCluster.FindSubClusterName(*point, level);
	vector<TString>::iterator itr = std::find( stringList.begin(), stringList.end(), nam);
	if ( itr == stringList.end()){
	  stringList.push_back(nam);
	}
	itr = std::find( stringList.begin(), stringList.end(), nam);
	hTTH->Fill(itr-stringList.begin(), point->fW);
  }
  vector<TH1*> bkgs;
  bkgs.push_back(hTTbar);
  bkgs.push_back(hTTW  );  
  MakeSimpleCard card(hTTH, bkgs, Form("datacard_Level_%d",level), 37000., false);
  card.doCard();
  gSystem->Exec(Form("combine -M Asymptotic -m %d datacard_Level_%d.txt",level,level));
  TFile* output = TFile::Open(Form("higgsCombineTest.Asymptotic.mH%d.root",level));
  TTree* Tlimit = (TTree*) output->Get("limit");
  Double_t limit = 0;
  Float_t quantileExpected = 0;
  Tlimit->SetBranchAddress("limit",&limit);
  Tlimit->SetBranchAddress("quantileExpected", &quantileExpected);
  for (Int_t entry = 0; entry < Tlimit->GetEntries(); entry++){
    Tlimit->GetEntry(entry);
    if (quantileExpected > 0.49 && quantileExpected < 0.51){
      cout << "Expected upper limit for level " << level <<  "  " << limit << endl;
      break;
    }
  }
    
  return limit;
}


void TargettedClustering::SignificancesEachLevel()
{
  TH1F* h = new TH1F("h","",7, 0.5,7.5);
  for (int k = 1; k < 8; ++k){
    h->SetBinContent(k, SignificanceAtLevel(k));
    h->SetBinError(k,0.);
  }
  
  h->SetMarkerStyle(kFullCircle);
  h->Draw("P");

}
