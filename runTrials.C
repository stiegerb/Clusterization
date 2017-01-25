void runTrials(int trial){

  // Run me with:
  // root -l -b runTrials.C

  gROOT->ProcessLine(".L MakeSimpleCard.C+");
  gROOT->ProcessLine(".L Significance.C+");
  gROOT->ProcessLine(".L TargettedClustering.C+");


  gROOT->ProcessLine("TargettedClustering* g = new TargettedClustering(2,2,0)");
  
  //for(int trial=1; trial<1001; ++trial)
  //  {
      gROOT->ProcessLine("delete g");
      gROOT->ProcessLine(Form("g = new TargettedClustering(2,2,%d)", trial));
      gROOT->ProcessLine("g->Test()");
      gROOT->ProcessLine("g->SignificancesEachLevel()");
      //  }
  
}
