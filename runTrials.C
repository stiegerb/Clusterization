void runTrials(int trial, TString module){

  // Run me with:
  // root -l -b runTrials.C

  cout << "Running trial " << trial << " for module " << module << endl;

  gROOT->ProcessLine(".L MakeSimpleCard.C+");
  gROOT->ProcessLine(".L Significance.C+");
  gROOT->ProcessLine(Form(".L %s.C+",module.Data()));
  
  
  gROOT->ProcessLine(Form("%s* g",module.Data()));
  
  //for(int trial=1; trial<1001; ++trial)
  //  {
  gROOT->ProcessLine(Form("g = new %s(2,3,%d)",module.Data(),trial));
  gROOT->ProcessLine("g->Test()");
  if(module=="TargettedClustering") gROOT->ProcessLine("g->SignificancesEachLevel()");
  gROOT->ProcessLine("g->VoronoiPlot()");
  gROOT->ProcessLine("delete g");
  //  }
  gApplication->Terminate(0);
}
