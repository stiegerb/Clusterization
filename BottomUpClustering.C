#include "BottomUpClustering.h"
#include "MakeSimpleCard.h"
#include "tdrstyle.C"
#include <TFile.h>
#include "TROOT.h"
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
#include <TLegend.h>
#include <TLatex.h>
#include <stdlib.h>


BottomUpClustering::BottomUpClustering(Int_t nLep, TString fileType, Int_t trial):
  nLep_(nLep),
  trial_(trial),
  fileType(fileType)
{
  ReadFromFiles();
  StartTheThing();

}

void BottomUpClustering::ReadFromFiles()
{
  if ( fileType == "txt"){
    readFromTxTFiles();
  }
  else if (fileType == "root"){
    readFromRootFiles();
  }
  else{
    cout << "Please select a correct file type" << endl;
    cout << "Script will crash in..." << endl;
    cout << "...3" << endl;
    sleep(1);
    cout << "...2" << endl;
    sleep(1);
    cout << "...1" << endl;
    sleep(1);
  }
}

void BottomUpClustering::StartTheThing()
{
  MakeFineBinning();
}


void BottomUpClustering::readFromRootFiles()
{
  // Todo:
  // Change file and branch names
  fTTbar.clear(); fTTH.clear(); fTTW.clear();
  TFile* f    = 0;
  TTree* tree = 0;
  TBranch *bX = 0;
  TBranch *bY = 0;
  TBranch *bW = 0;

  Double_t x = 0.;
  Double_t y = 0.;
  Double_t w = 0.;

  // Reading ttbar
  f = TFile::Open( (nLep_ == 3) ? "data/ttbar3l.root" : "data/ev_2lss_TT.root");
  tree = (TTree*) f->Get("tree");
  tree->Branch("kinMVA_2lss_ttbar_withBDTv8", &x, "kinMVA_2lss_ttbar_withBDTv8/F");
  tree->Branch("kinMVA_2lss_ttbar_withBDTv8", &y, "kinMVA_2lss_ttV_withHj/F");
  tree->Branch("_weight_", &w, "_weight_/F");
  for (int entr = 0; entr < tree->GetEntries(); entr++){
    tree->GetEntry(entr);
    Point point = Point(x,y,2*w);
    if (entr % 2 == 0)
      fTTbar.push_back(point);
    else
      fTTbarMC.push_back(point);
  }
  f->Close();
  cout << "TTbar events " << fTTbar.size() << endl;


  // Reading ttH
  f = TFile::Open( (nLep_ == 3) ? "data/tth3l.root" : "data/ev_2lss_TTHnobb_pow.root");
  tree = (TTree*) f->Get("tree");
  tree->Branch("kinMVA_2lss_ttbar_withBDTv8", &x, "kinMVA_2lss_ttbar_withBDTv8/F");
  tree->Branch("kinMVA_2lss_ttbar_withBDTv8", &y, "kinMVA_2lss_ttV_withHj/F");
  tree->Branch("_weight_", &w, "_weight_/F");
  for (int entr = 0; entr < tree->GetEntries(); entr++){
    tree->GetEntry(entr);
    Point point = Point(x,y,2*w);
    if (entr % 2 == 0)
      fTTH.push_back(point);
    else
      fTTHMC.push_back(point);
  }
  f->Close();
  cout << "TTH events " << fTTH.size() << endl;

  // Reading ttW
  f = TFile::Open( (nLep_ == 3) ? "data/ttw3l.root" : "data/ev_2lss_TTV.root");
  tree = (TTree*) f->Get("tree");
  tree->Branch("kinMVA_2lss_ttbar_withBDTv8", &x, "kinMVA_2lss_ttbar_withBDTv8/F");
  tree->Branch("kinMVA_2lss_ttbar_withBDTv8", &y, "kinMVA_2lss_ttV_withHj/F");
  tree->Branch("_weight_", &w, "_weight_/F");
  for (int entr = 0; entr < tree->GetEntries(); entr++){
    tree->GetEntry(entr);
    Point point = Point(x,y,2*w);
    if (entr % 2 == 0)
      fTTW.push_back(point);
    else
      fTTWMC.push_back(point);
  }
  f->Close();
  cout << "TTW events " << fTTW.size() << endl;

}

void BottomUpClustering::readFromTxTFiles()
{


  fTTbar.clear(); fTTH.clear(); fTTW.clear();
  ifstream f;
  nLep_==3 ? f.open("data/ttbar3l.txt") : f.open("data/ttbar.txt");
  Int_t count = 0;
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> y >> x >> w; //
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

  nLep_==3 ? f.open("data/tth3l.txt") : f.open("data/tth.txt");
  while (f){
    Double_t x(0), y(0), w(0), z(0);
    f >> y >> x >> w >> z;
    if (w == 0) continue;
    if ( TMath::Abs(x) > 2.) continue;
    if ( TMath::Abs(y) > 2.) continue;

    w *= z;
    Point point = Point(x,y,2*w);
    if (count%2 == 0)
      fTTH.push_back(point);
    else
      fTTHMC.push_back(point);
    count++;
  }
  f.close();
  cout << "TTH events " << fTTH.size() << endl;

  nLep_==3 ? f.open("data/ttv3l.txt") : f.open("data/ttv.txt");
  while (f){
    Double_t x = 0; Double_t y = 0; Double_t w = 0;
    f >> y >> x >> w;
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
  cout << "TTW events " << fTTW.size() << endl;
}

void BottomUpClustering::MakeFineBinning()
{
  hFineBinning = new TH2F("hFineBinning","",20,-1.,1.,20,-1.,1.);
  int counter = 0;
  for (int i = 0; i < hFineBinning->GetXaxis()->GetNbins(); ++i){
    for (int j = 0; j < hFineBinning->GetYaxis()->GetNbins(); ++j){
      hFineBinning->SetBinContent(hFineBinning->GetBin(i+1,j+1),counter);
      vector<Int_t> vct; vct.push_back(hFineBinning->GetBin(i+1,j+1));
      clusters.push_back(vct);
      counter++;
    }
  }
  // we start from the fine binning
  hTargetBinning = (TH2F*) hFineBinning->Clone("hTargetBinning");
  return;
}

void BottomUpClustering::ReMakeTarget()
{
  for (unsigned int k = 0; k < clusters.size(); ++k){
    for (auto& minibin : clusters[k]){ // loop over its minibins (minibin is the bin number)
      hTargetBinning->SetBinContent( minibin, k);
    }
  }  
  
  return;
}


void BottomUpClustering::GetEventsInCluster(Int_t k, Double_t& ttbar, Double_t& ttv, Double_t& tth,Double_t& ttbar_e, Double_t& ttv_e, Double_t& tth_e)
{
  TH1F* hTTbar = new TH1F("hTTbar","",clusters.size(),-0.5,clusters.size()-0.5);
  TH1F* hTTW   = new TH1F("hTTW"  ,"",clusters.size(),-0.5,clusters.size()-0.5);
  TH1F* hTTH   = new TH1F("hTTH"  ,"",clusters.size(),-0.5,clusters.size()-0.5);
  
  for (auto& pt : fTTbar){
    hTTbar->Fill(hTargetBinning->GetBinContent( hTargetBinning->FindBin(pt.fX, pt.fY) ), pt.fW);
  }
  for (auto& pt : fTTW){
    hTTW->Fill(hTargetBinning->GetBinContent( hTargetBinning->FindBin(pt.fX, pt.fY) ), pt.fW);
  }
  for (auto& pt : fTTH){
    hTTH->Fill(hTargetBinning->GetBinContent( hTargetBinning->FindBin(pt.fX, pt.fY) ), pt.fW);
  }

  ttbar = hTTbar -> GetBinContent( hTTbar -> FindBin(k) );
  ttv   = hTTW   -> GetBinContent( hTTW   -> FindBin(k) );
  tth   = hTTH   -> GetBinContent( hTTH   -> FindBin(k) );
  ttbar_e = hTTbar -> GetBinError( hTTbar -> FindBin(k) );
  ttv_e   = hTTW   -> GetBinError( hTTW   -> FindBin(k) );
  tth_e   = hTTH   -> GetBinError( hTTH   -> FindBin(k) );
		     
  delete hTTbar;
  delete hTTW;
  delete hTTH;
}

Double_t BottomUpClustering::GetFOM(Int_t k1, Int_t k2)
{
  Double_t ttbar1, ttv1, tth1;
  Double_t ttbar2, ttv2, tth2;
  Double_t ttbar1_e, ttv1_e, tth1_e;
  Double_t ttbar2_e, ttv2_e, tth2_e;
  GetEventsInCluster(k1,ttbar1, ttv1, tth1,ttbar1_e, ttv1_e, tth1_e);
  GetEventsInCluster(k2,ttbar2, ttv2, tth2,ttbar2_e, ttv2_e, tth2_e);
  Double_t propTTbar1 = ttbar1/(ttbar1 + ttv1 + tth1);
  Double_t propTTbar2 = ttbar2/(ttbar2 + ttv2 + tth2);
  Double_t propTTV1   = ttv1/(ttbar1 + ttv1 + tth1);
  Double_t propTTV2   = ttv2/(ttbar2 + ttv2 + tth2);

  // cout << ttbar1 << " " << ttv1 << " " << tth1 << " " <<  ttbar2 << " " << ttv2 << " " << tth2 <<endl;
  Double_t propTTbar1_e = ttbar1_e / (ttbar1 + ttv1 + tth1) + ttbar1*(ttv1_e+tth1_e)/((ttbar1 + ttv1 + tth1)*(ttbar1 + ttv1 + tth1));
  Double_t propTTbar2_e = ttbar2_e / (ttbar2 + ttv2 + tth2) + ttbar2*(ttv2_e+tth2_e)/((ttbar2 + ttv2 + tth2)*(ttbar2 + ttv2 + tth2));
  Double_t propTTV1_e   = ttv1_e   / (ttbar1 + ttv1 + tth1) + ttv1*(ttbar1_e+tth1_e)/((ttbar1 + ttv1 + tth1)*(ttbar1 + ttv1 + tth1));
  Double_t propTTV2_e   = ttv2_e   / (ttbar2 + ttv2 + tth2) + ttv2*(ttbar2_e+tth2_e)/((ttbar2 + ttv2 + tth2)*(ttbar2 + ttv2 + tth2));
  // if the fom is close to zero, the clusters are the same
  if ( (ttbar1*ttv1*tth1*ttbar2*ttv2*tth2) == 0) return 0.;
  return TMath::Power( (propTTbar1 - propTTbar2) / (propTTbar1_e + propTTbar2_e), 2) + TMath::Power( (propTTV1 - propTTV2) / (propTTV1_e + propTTV2_e), 2) + 37000*2.*(ttbar1+ ttv1+ tth1+ttbar2+ ttv2+ tth2);

}

void BottomUpClustering::Recluster()
{
  for (int it = 0; it < 390; ++it){
    Double_t minFOM = 1000;
    std::pair<Int_t,Int_t> bestPair; 
    int index1 = 0;
    int index2 = 0;
    cout << "We have " << clusters.size() << " clusters" << endl;

    for (unsigned int k1 = 0; k1 < clusters.size(); ++k1){
      for (unsigned int k2 = 0; k2 < clusters.size(); ++k2){
	if (k2 <= k1) continue;
	if (!AreClustersTogether(k1,k2)) continue;
	Double_t fom = GetFOM(k1,k2);
	if (fom < minFOM){
	  bestPair = std::make_pair(k1,k2);
	  minFOM = fom;
	  if (fom == 0) break;
	}
      }
    }
    cout << "Best reclustering pair is " << bestPair.first << " " << bestPair.second <<  " with fom " << minFOM << endl;
    if (bestPair.first == bestPair.second) break;
    vector<vector<Int_t>> newClusters;
    for (unsigned int k = 0; k< clusters.size(); ++k){
      if (k == bestPair.second){
	// cout << bestPair.first << " " << newClusters.size() << endl;
	// cout << clusters.size() << endl;
	newClusters[bestPair.first].insert(newClusters[bestPair.first].end(),clusters[k].begin(),clusters[k].end());
      }
      else{
	newClusters.push_back(clusters[k]);
      }
    }
    clusters = newClusters;
    ReMakeTarget();
  }
  TCanvas* c1 = new TCanvas();
  hTargetBinning->Draw("colz text");
  TCanvas* c2 = new TCanvas();
  Test();
  return;
}

Double_t BottomUpClustering::GetCluster(Point pt)
{
  return hTargetBinning->GetBinContent( hTargetBinning->FindBin(pt.fX, pt.fY));
}

void BottomUpClustering::Test()
{
  TCanvas* c = new TCanvas();
  c->cd();
  setTDRStyle();

  TH1F* hTTbar = new TH1F("hTTbar","",clusters.size(), -0.5, clusters.size()-0.5);
  TH1F* hTTW   = new TH1F("hTTW"  ,"",clusters.size(), -0.5, clusters.size()-0.5);
  TH1F* hTTH   = new TH1F("hTTH"  ,"",clusters.size(), -0.5, clusters.size()-0.5);
  THStack* mc  = new THStack("mc","mc");
  vector<Point>::iterator point;
  cout << "TTbar size is " << fTTbarMC.size() << endl;
  for (point = fTTbarMC.begin(); point != fTTbarMC.end(); ++point)
    hTTbar->Fill( GetCluster( *point), 36500*point->fW);
  for (point = fTTWMC.begin(); point != fTTWMC.end(); ++point)
    hTTW->Fill( GetCluster( *point), 36500*point->fW);
  for (point = fTTHMC.begin(); point != fTTHMC.end(); ++point)
    hTTH->Fill( GetCluster( *point), 36500*point->fW);
  cout << hTTbar->Integral() << " " << hTTbar->GetEntries() << endl;
  hTTbar->SetFillColor( kRed     );
  hTTH->SetFillColor( kBlue    );
  hTTW->SetFillColor( kMagenta );

  mc->Add( hTTbar ); mc->Add(hTTW); mc->Add(hTTH);
  mc->Draw("HIST");
  mc->SetMaximum(1.5* mc->GetMaximum());
  mc->GetHistogram()->GetYaxis()->SetTitle("Expected events/bin");
  mc->GetHistogram()->GetXaxis()->SetTitle("Bin in the bdt1#times bdt2 plane");
  mc->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
  mc->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);
  mc->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
  mc->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);

  TLegend* l = new TLegend(0.7,0.8,0.9,0.9);
  l->AddEntry(hTTH  , "ttH signal", "f");
  l->AddEntry(hTTW  , "ttV"       , "f");
  l->AddEntry(hTTbar, "tt"        , "f");
  l->Draw();

  TLatex latex;
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  latex.SetTextFont(62);
  latex.DrawLatexNDC(.1,.95,"CMS Simulation");
  latex.DrawLatexNDC(.7,.95,"#it{36 fb^{-1}}");

  c->Modified();
  c->Update();


  // c->Print(Form("recursiveNoOrdering_%dl_trial%d.png",nLep_,trial_));
  // c->Print(Form("recursiveNoOrdering_%dl_trial%d.pdf",nLep_,trial_));




  // vector<Double_t> ttbarLikeness;
  // vector<Double_t> ttwLikeness;
  // for (int k = 0; k < gIndex; ++k){
  //   ttbarLikeness.push_back( hTTbar->GetBinContent(k+1) / (hTTbar->GetBinContent(k+1) + hTTW->GetBinContent(k+1) + hTTH->GetBinContent(k+1)));
  //   ttwLikeness.push_back( hTTW->GetBinContent(k+1) / (hTTbar->GetBinContent(k+1) + hTTW->GetBinContent(k+1) + hTTH->GetBinContent(k+1)));
  // }
  // TCanvas* c1 = new TCanvas();
  // TGraph* gr = new TGraph(gIndex, &ttbarLikeness[0], &ttwLikeness[0]);
  // gr->SetMarkerStyle(kFullCircle);
  // gr->Draw("A,P");


}

Int_t BottomUpClustering::SortedThing(Int_t bin)
{
  for (unsigned int k = 0; k < SoverB.size(); ++k){
    if ( SoverB[k].second == bin ) return k;
  }
  cout << "[SortedThing::WTF]" << endl;
  return -1;
}

void BottomUpClustering::StoreToFile()
{
  TFile* binning = TFile::Open(Form("binning_%dl.root",nLep_),"recreate");
  hTargetBinning->Write();
  binning->Close();
}

void BottomUpClustering::LoadClusterFromFile()
{
  TFile* binning = TFile::Open(Form("binning_%dl.root",nLep_),"recreate");
  hTargetBinning = (TH2F*) binning->Get("hTargetBinning");
  hTargetBinning->SetDirectory(0);
  binning->Close();

  Int_t nClusters = hTargetBinning->GetMaximum();
  for (int k = 0; k < nClusters; ++k){
    vector<Int_t> kk; // empty cluster
    clusters.push_back(kk);
  }
  for (int i = 1; i < hTargetBinning->GetXaxis()->GetNbins()+1; ++i){
      for (int j = 1; j < hTargetBinning->GetXaxis()->GetNbins()+1; ++j){
	clusters[ int(hTargetBinning->GetBinContent( hTargetBinning->GetBin(i,j)))].push_back(hTargetBinning->GetBin(i,j));
      }
  }

}

//   // Final significance per bin
//   double indexes[fSignificances.size()];
//   double signifs[fSignificances.size()];
//   for(unsigned int idx=0; idx<fSignificances.size(); ++idx)
//     {
//       indexes[idx]=idx;
//       signifs[idx]=fSignificances[idx];
//     }
//   TGraph* g = new TGraph(fSignificances.size(), indexes, signifs);
//   g->SetName("significancePerBin_finalSet");
//   g->Write();

//   binning->Close();
// }


// void BottomUpClustering::VoronoiPlot()
// {
//   vector<TGraph*> graphs; graphs.clear();
//   vector<Double_t>* X = new vector<Double_t>[clusters.size()];
//   vector<Double_t>* Y = new vector<Double_t>[clusters.size()];

//   cout << "Calculating points" << endl;
//   for (Double_t x = -1; x < 1.; x = x + 1e-3){

//       for (Double_t y = -1; y < 1.; y = y + 1e-3){
// 	Int_t k = 0; //mainCluster.FindUnclusterizableCluster(Point(x,y,-1));
// 	X[k].push_back(x);
// 	Y[k].push_back(y);
//       }
//   }



//   TCanvas* c = new TCanvas();
//   c->cd();
//   setTDRStyle();

//   TH1F* hDummy = new TH1F("hDummy","",2,-1,1);
//   hDummy->SetBinContent(1, 1.);
//   hDummy->SetBinContent(2,-1.);
//   hDummy->SetLineColor(kWhite);
//   hDummy->GetYaxis()->SetRangeUser(-1.,1.);
//   hDummy->GetXaxis()->SetTitle("BDT(ttH,tt)");
//   hDummy->GetYaxis()->SetTitle("BDT(ttH,ttV)");
//   hDummy->Draw();
//   cout << "Done... now plotting" << endl;

//   TText t;
//   t.SetTextSize(0.08);
//   for (unsigned int k = 0; k < clusters.size(); ++k){
//     graphs.push_back(new TGraph( X[k].size(), &X[k][0], &Y[k][0] ));
//     graphs[k]->SetMarkerColor(k);
//     graphs[k]->SetMarkerStyle(6);
//     graphs[k]->Draw("PSAME");

//     t.SetTextColor(k+1);
//     //    t.DrawText(fCentroids[k].fX, fCentroids[k].fY, Form("%d",k));

//   }


//   c->Print(Form("voronoi_%dl_trial%d.png",nLep_,trial_));
//   if(trial_==0) // Only save PDF when running in single mode. Otherwise, huuuuge clogging of backup server
//     c->Print(Form("voronoi_%dl_trial%d.pdf",nLep_,trial_));
//   cout << "REACTIVATE VORONOI.PDF" << endl;
// }




// void BottomUpClustering::SortedTest()
// {

//   TCanvas* c = new TCanvas();
//   c->cd();
//   setTDRStyle();

//   TH1F* hTTbar = new TH1F("hTTbar","",clusters.size(), -0.5, clusters.size()-0.5);
//   TH1F* hTTW   = new TH1F("hTTW"  ,"",clusters.size(), -0.5, clusters.size()-0.5);
//   TH1F* hTTH   = new TH1F("hTTH"  ,"",clusters.size(), -0.5, clusters.size()-0.5);


//   THStack* mc  = new THStack("mc","mc");
//   vector<Point>::iterator point;

//   for (point = fTTbarMC.begin(); point != fTTbarMC.end(); ++point)
//     {
//       hTTbar->Fill( SortedThing( mainCluster.FindUnclusterizableCluster( *point)), 36500*point->fW);
//     }
//   for (point = fTTWMC.begin(); point != fTTWMC.end(); ++point)
//     {
//       hTTW->Fill( SortedThing( mainCluster.FindUnclusterizableCluster( *point)), 36500*point->fW);
//     }
//   for (point = fTTHMC.begin(); point != fTTHMC.end(); ++point)
//     {
//       hTTH->Fill( SortedThing( mainCluster.FindUnclusterizableCluster( *point)), 36500*point->fW);
//     }
//   hTTbar->SetFillColor( kRed     );
//   hTTH->SetFillColor( kBlue    );
//   hTTW->SetFillColor( kMagenta );


//   cout << "Now simple card test " << endl;
//   vector<TH1*> bkgs;
//   bkgs.push_back(hTTbar);
//   bkgs.push_back(hTTW  );

//   MakeSimpleCard card(hTTH, bkgs, "datacard_recursiveclustering", 1., false);
//   card.doCard();




//   mc->Add( hTTbar ); mc->Add(hTTW); mc->Add(hTTH);
//   mc->Draw("HIST");

//   mc->GetHistogram()->GetYaxis()->SetTitle("Expected events/bin");
//   mc->GetHistogram()->GetXaxis()->SetTitle("Bin in the bdt1#times bdt2 plane");
//   mc->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
//   mc->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);
//   mc->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
//   mc->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);

//   TLegend* l = new TLegend(0.7,0.8,0.9,0.9);
//   l->AddEntry(hTTH  , "ttH signal", "f");
//   l->AddEntry(hTTW  , "ttV"       , "f");
//   l->AddEntry(hTTbar, "tt"        , "f");
//   l->Draw();

//   TLatex latex;
//   latex.SetTextSize(0.05);
//   latex.SetTextAlign(13);  //align at top
//   latex.SetTextFont(62);
//   latex.DrawLatexNDC(.1,.95,"CMS Simulation");
//   latex.DrawLatexNDC(.7,.95,"#it{36 fb^{-1}}");

//   c->Modified();
//   c->Update();


//   c->Print(Form("recursive_%dl.png",nLep_));
//   c->Print(Form("recursive_%dl.pdf",nLep_));



// }




bool BottomUpClustering::AreClustersTogether(int k1, int k2)
{
  // cout << "Checking if clusters " << k1 << " " << k2 << " are together" << endl;
  // cout << "Size " << clusters[k1].size() << endl;
  for (auto& bin : clusters[k1]){
    Int_t x,y,z;
    hTargetBinning->GetBinXYZ(bin, x,y,z);
    Double_t x_  = hTargetBinning->GetXaxis()->GetBinCenter(x);
    Double_t x_d = hTargetBinning->GetXaxis()->GetBinWidth(x) ;
    Double_t y_  = hTargetBinning->GetYaxis()->GetBinCenter(y);
    Double_t y_d = hTargetBinning->GetYaxis()->GetBinWidth(y) ;

    if ( x_+x_d < 1){
      if ( hTargetBinning->GetBinContent(hTargetBinning->FindBin(x_+x_d,y_))==k2){
	return true;
      }
    }
    if ( x_-x_d > -1){
      if ( hTargetBinning->GetBinContent(hTargetBinning->FindBin(x_-x_d,y_))==k2){
	return true;
      }
    }
    if (y_+y_d < 1){
      if ( hTargetBinning->GetBinContent(hTargetBinning->FindBin(x_,y_+y_d))==k2){
	return true;
      }
    }
    if (y_-y_d > -1){
      if ( hTargetBinning->GetBinContent(hTargetBinning->FindBin(x_,y_-y_d))==k2){
	return true;
      }
    }
  }

  return false;
  // for (auto& gr : contours){
  //   gr->SetLineColor(kRed);
  //   gr->SetLineWidth(2);
  //   gr->Draw("LSAME");
  //   cout << "Checking new contour" << endl;
  //   bool isK1Limit = false;
  //   bool isK2Limit = false;
  //   for (int point = 1; point < gr->GetN()+1; ++point){
  //     gr->GetPoint(point,x,y);
  //     //cout << x << " " << y <<  " " << hTargetBinning->GetBinContent(hTargetBinning->FindBin(x,y)) << endl;
  //     for (float phi = 0; phi < 2*TMath::Pi(); phi = phi+0.001){
  // 	Point point(x+0.001*TMath::Cos(phi),y+0.001*TMath::Sin(phi),-1);
  // 	bool isK1 = (TMath::Abs(hTargetBinning->GetBinContent( hTargetBinning->FindBin(x,y)) - k1) < 0.01);
  // 	bool isK2 = (TMath::Abs(hTargetBinning->GetBinContent( hTargetBinning->FindBin(x,y)) - k2) < 0.01);
  // 	// if (isK1 && !isK1Limit) cout << "Point " << point << "is in cluster " << k1 << endl;
  // 	// if (isK2 && !isK2Limit) cout << "Point " << point << "is in cluster " << k2 << endl;
  // 	isK1Limit = isK1Limit || isK1;
  // 	isK2Limit = isK2Limit || isK2;
  //     }
  //   }
  //   if (isK1Limit && isK2Limit) return true;
  // }
  return false; 
}



void BottomUpClustering::TestBoundaries()
{

  hTargetBinning->Draw("colz text");
  for (unsigned int j = 0; j < clusters.size(); ++j){
    for (unsigned int k = 0; k < clusters.size(); ++k){
      if (AreClustersTogether(j,k)) cout << "They are together "  << j << " " << k << endl;
    }
  }

  return;
}

Int_t BottomUpClustering::NearestClusters(Int_t bin)
{
  Int_t x,y,z;
  hTargetBinning->GetBinXYZ(bin, x,y,z);
  Double_t x_  = hTargetBinning->GetXaxis()->GetBinCenter(x);
  Double_t y_  = hTargetBinning->GetYaxis()->GetBinCenter(y);
  Double_t dist = 999;
  Int_t theCluster = -1;
  for (unsigned int cl = 0; cl < clusters.size(); ++cl){
    for (auto& b: clusters[cl]){
      hTargetBinning->GetBinXYZ(b, x,y,z);
      Double_t x2  = hTargetBinning->GetXaxis()->GetBinCenter(x);
      Double_t y2  = hTargetBinning->GetYaxis()->GetBinCenter(y);
      Double_t d   = (x_-x2)*(x_-x2)+(y_-y2)*(y_-y2);
      if (d < dist){
	dist   = d;
	theCluster = cl;
      }
    }
  }
  if (theCluster < 0){
    cout << "Something wrong is happening in nearestclusters" << endl;
    cout << "Will return -1 and it will crash" << endl;
  }
  return theCluster;

}

void BottomUpClustering::ReCleanClusters()
{
  for (auto& bin : clusters[0]){
    clusters[NearestClusters(bin)].push_back(bin);
  }
  clusters.erase(clusters.begin());
  ReMakeTarget();
  Test();

}
