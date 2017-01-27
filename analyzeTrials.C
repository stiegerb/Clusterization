#include <iostream>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2D.h>
#include <TStyle.h>

using namespace std;

void analyzeTrials(int nTrials){

  gStyle->SetOptStat(0);

  std::vector<TGraph*> g; g.clear();
  std::vector<TGraph*> gnbins; gnbins.clear();
  TFile* output=nullptr;
  TFile* shapes=nullptr;
  TTree* Tlimit=nullptr;
  TH1* nbins=nullptr;
  TH2D* dep = new TH2D("dep", "Number of trials per bin;Level;Nbins", 7, 0.5, 7.5, 20, 0.5, 20.5);
  double levels[7] = {1., 2., 3., 4., 5., 6., 7.};

  double minLimit(99999.);
  double minLabel(-1);
  for(int trial=1; trial<nTrials; ++trial)
    {
      cout << "[" << trial << "/" << nTrials << "]" << endl;
      double limits[7];
      double nb[7];
      for(int level=1; level<8; ++level)
        {
          output = TFile::Open(Form("higgsCombineTrial%d.Asymptotic.mH%d.root",trial,level));
          shapes = TFile::Open(Form("datacard%d_Level_%d_shapes.root",trial,level));
          
          if(output)
            {
              Tlimit = (TTree*) output->Get("limit");
              if(Tlimit)
                {
                  Double_t limit(0);
                  Tlimit->SetBranchAddress("limit",&limit);
                  Tlimit->GetEntry(2);
                  limits[level]=limit;
                  if(level==7)
                    {
                      if(limit<minLimit)
                        {
                          minLimit=limit;
                          minLabel=trial;
                        }
                    }
                  delete Tlimit;
                  output->Close();
                  delete output;
                }
            }
          if(shapes)
            {
              nbins = (TH1*) shapes->Get("data_obs");
              nb[level] = nbins->GetNbinsX();
              delete nbins;
              shapes->Close();
              delete shapes;
            }
        }
      g.push_back(new TGraph(7, levels, limits));
      g[trial-1]->SetName(Form("g%d",trial));
      gnbins.push_back(new TGraph(7, levels, nb));
      gnbins[trial-1]->SetName(Form("nbi%d",trial));
      for(size_t lev=0; lev<7; ++lev)
        {
          dep->Fill(levels[lev], nb[lev]);
        }
    }
  
  TCanvas* c = new TCanvas("c", Form("%d bootstrap trials", nTrials), 1000, 1000);
  c->cd();
  bool painted(false);
  int utrial(0);
  for(auto& graph : g)
    {
      if(!painted)
        {
          graph->Draw("APL");
          painted=true;
        }
      else
        graph->Draw("PL");
      cout << "Fetched graph for trial " << utrial << endl;
      utrial++;
    }
  c->Modified();
  c->Update();
  c->Print("expectedLimitTrials.png");
  c->Print("expectedLimitTrials.pdf");

  c->Clear();
  c->Update();
  painted=false;
  for(auto& graph : gnbins)
    {
      if(!painted)
        {
          graph->Draw("APL");
          graph->GetHistogram()->SetMaximum(20.);
          painted=true;
        }
      else
        graph->Draw("PL");
    }
  c->Modified();
  c->Update();
  c->Print("nbinsPerLevelTrials.png");
  c->Print("nbinsPerLevelTrials.pdf");

  c->Clear();
  c->Update();
  dep->Draw("contz");
  c->Modified();
  c->Update();
  c->Print("fulldependency.png");
  c->Print("fulldependency.pdf");


  cout << "Minimum is at " << minLimit << ", corresponding to trial " << minLabel << endl;
}
