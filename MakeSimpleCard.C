#include "MakeSimpleCard.h"



MakeSimpleCard::MakeSimpleCard(TH1* sig, vector<TH1*> bkg, TString cardName, double lumi, bool debug):
  sig_(sig),
  bkg_(bkg),
  cardName_(cardName),
  lumi_(lumi),
  systSig_(0.5),
  systBkg_(0.10),
  nBkg_(bkg.size()+1),
  debug_(debug)
{
  cout << "[MakeSimpleCard::MakeSimpleCard] Welcome to MakeSimpleCard." << endl;
}

void MakeSimpleCard::doCard()
{
  cout << "[MakeSimpleCard::doCard] I will now prepare a simple card." << endl;
  card_.open(cardName_+".txt");
  cout << "Blsdjgshldh" << endl;
  card_ << "# This is a sample card for categorization study.\n";

  Renormalize();
  
  FillFirstBlock();
  
  FillObserved();

  FillRates();

  cout << "[MakeSimpleCard::doCard] Need to implement statistical uncertainties" << endl;

  WriteCard();

  cout << "[MakeSimpleCard::doCard] Please have a nice day." << endl;
}

void MakeSimpleCard::Renormalize()
{
  sig_->Scale(lumi_);
  for(auto& bkg : bkg_)
    bkg->Scale(lumi_);
  if(debug_) cout << "[MakeSimpleCard::Renormalize] Method finished successfully." << endl;
}

void MakeSimpleCard::FillFirstBlock()
{
  card_ << "imax 1\n";
  card_ << "jmax 1\n";
  card_ << "kmax *\n";
  card_ << "---------------\n";
  card_ << "shapes * * " << cardName_ << "_shapes.root $PROCESS $PROCESS_$SYSTEMATIC\n";
  card_ << "---------------\n";
  if(debug_) cout << "[MakeSimpleCard::FillFirstBlock] Method finished successfully." << endl;
}

void MakeSimpleCard::FillObserved()
{
  data_=nullptr;
  double datayield(0.);
  for(auto& bkg : bkg_)
    {
      if(debug_) cout << "[MakeSimpleCard::FillObserved] Looping. " << bkg->GetName() << " " << bkg->Integral() << endl;
      datayield += bkg->Integral();
      if(data_) data_->Add(bkg,1);
      else      data_= (TH1*) bkg->Clone("data_obs");
    }
  if(debug_) cout << "[MakeSimpleCard::FillObserved] Ended loop." << endl;
  card_ << "bin a\n";
  card_ << "observation " << data_->Integral() << "\n";
  card_ << "------------------------------\n";
  if(debug_) cout << "[MakeSimpleCard::FillObserved] Method finished successfully." << endl;
}

void MakeSimpleCard::FillRates()
{
  card_ << "bin             a          ";
  for(auto& bkg : bkg_)
    card_ << "      a";
  card_<< "\n";

  card_ << "process         signal     ";
  for(auto& bkg : bkg_)
    card_ << "       " << bkg->GetName();
  card_ << "\n";
  
  size_t iSample(0);
  card_ << "process         " << iSample << "          ";
  for(auto& bkg : bkg_)
    {
      iSample++;
      card_ << "       " << iSample;
    }
  card_ << "\n";

  card_ << "lumi        " << 0.02;
  for(auto& bkg : bkg_)
    {
      card_ << "       " << 0.02;
    }
  card_ << "\n";

  card_ << "sigXsec        " << double(1+systSig_);
  for(auto& bkg : bkg_)
    {
      card_ << "      - ";
    }
  card_ << "\n";

  for(auto& bkg : bkg_)
    {
      card_ << bkg->GetName() << "Xsec        - ";
      for(auto& ibkg : bkg_)
        {
          ibkg == bkg ? card_ << "      " << double(1+systBkg_) : card_ << "      -";
        }
      card_ << "\n";
    }
  card_ << "--------------------------------\n";
  if(debug_) cout << "[MakeSimpleCard::FillRates] Method finished successfully." << endl;
}

void MakeSimpleCard::WriteShapes()
{
  shapesFile_ = TFile::Open(cardName_+"_shapes.root", "RECREATE");
  
  sig_->Write();
  for(auto& bkg : bkg_)
    bkg->Write();
  
  shapesFile_->Close();
  if(debug_) cout << "[MakeSimpleCard::WriteShapes] Method finished successfully." << endl;
}

void MakeSimpleCard::WriteCard()
{
  card_.close();
  cout << "[MakeSimpleCard::writeCard] Card has been written." << endl;
  cout << "                            Process it in CMSSW_7_1_X with:" << endl;
  cout << "                            combine -M Asymptotic " << cardName_ << ".txt" << endl;
  if(debug_) cout << "[MakeSimpleCard::WriteCard] Method finished successfully." << endl;
}
