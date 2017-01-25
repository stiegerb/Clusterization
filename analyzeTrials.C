#include <iostream>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TCanvas.h>

using namespace std;

void analyzeTrials(int nTrials){

  std::vector<TGraph*> g; g.clear();
  TFile* output=nullptr;
  TTree* Tlimit=nullptr;
  double levels[7] = {1., 2., 3., 4., 5., 6., 7.};

  double minLimit(99999.);
  double minLabel(-1);
  for(int trial=1; trial<nTrials; ++trial)
    {
      double limits[7];
      for(int level=1; level<8; ++level)
        {
          output = TFile::Open(Form("higgsCombineTrial%d.Asymptotic.mH%d.root",trial,level));
          if(!output) continue;
          Tlimit = (TTree*) output->Get("limit");
          if(!Tlimit) continue;
          Double_t limit(0);
          Tlimit->SetBranchAddress("limit",&limit);
          Tlimit->GetEntry(3);
          limits[level]=limit;
          cout << "Level: " << level << ", limit: " << limit << endl;
          if(level==7)
            {
              if(limit<minLimit)
                {
                  minLimit=limit;
                  minLabel=trial;
                }
            }
          output->Close();
          delete output;
        }
      g.push_back(new TGraph(7, levels, limits));
      cout << "Fetched value for trial " << trial << endl;
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
  c->Print("bootstrapTrials.png");
  c->Print("bootstrapTrials.pdf");


  cout << "Minimum is at " << minLimit << ", corresponding to trial " << minLabel << endl;
}
