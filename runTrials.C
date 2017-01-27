void runTrials(int trial){

  // Run me with:
  // root -l -b runTrials.C

  cout << "Running trial " << trial << endl;

  gROOT->ProcessLine(".L MakeSimpleCard.C+");
  gROOT->ProcessLine(".L Significance.C+");
  gROOT->ProcessLine(".L TargettedClustering.C+");


  gROOT->ProcessLine("TargettedClustering* g");
  
  //for(int trial=1; trial<1001; ++trial)
  //  {
      gROOT->ProcessLine(Form("g = new TargettedClustering(2,3,%d)", trial));
      gROOT->ProcessLine("g->Test()");
      gROOT->ProcessLine("g->SignificancesEachLevel()");
      gROOT->ProcessLine("g->VoronoiPlot()");
      gROOT->ProcessLine("delete g");
      //  }
      gApplication->Terminate(0);
}
