#include "MakeSimpleCard.h"

#include <sstream>

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

  WriteShapes();

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
  card_ << "imax *\n";
  card_ << "jmax *\n";
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

  card_ << "process         " << sig_->GetName() << "      ";
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

  card_ << "rate        " << sig_->Integral();
  for(auto& bkg : bkg_)
    {
      card_ << "       " << bkg->Integral();
    }
  card_ << "\n";

  card_ << "-------------------\n";

  card_ << "lumi    lnN    " << 1.02;
  for(auto& bkg : bkg_)
    {
      card_ << "       " << 1.02;
    }
  card_ << "\n";

  card_ << "sigXsec lnN       " << double(1+systSig_);
  for(auto& bkg : bkg_)
    {
      card_ << "      - ";
    }
  card_ << "\n";

  for(auto& bkg : bkg_)
    {
      card_ << bkg->GetName() << "Xsec lnN       - ";
      for(auto& ibkg : bkg_)
        {
          ibkg == bkg ? card_ << "      " << double(1+systBkg_) : card_ << "      -";
        }
      card_ << "\n";
    }
  card_ << "--------------------------------\n";
  // Now fill up statistical uncertainties
  ostringstream convert;   // stream used for bin name conversion  
  for(int ibin=1; ibin<=sig_->GetNbinsX(); ++ibin)
    {
      convert << ibin;      // insert the textual representation of 'Number' in the characters in the stream
      TString binName(convert.str());
      card_ << "bin" << binName << sig_->GetName() << "Stat shape    1   ";
      for(auto& bkg : bkg_)
        {
          card_ << "      - ";
        }
      card_ << "\n";
      for(auto& bkg : bkg_)
        {
          card_ << "bin" << binName << bkg->GetName() << "Stat shape     -  ";
          for(auto& ibkg : bkg_)
            {
              ibkg == bkg ? card_ << "     1 " : card_ << "     - ";
            }
          card_ << "\n";
        }
    }
  card_ << "--------------------------------\n";
  if(debug_) cout << "[MakeSimpleCard::FillRates] Method finished successfully." << endl;
}

void MakeSimpleCard::WriteShapes()
{
  shapesFile_ = TFile::Open(cardName_+"_shapes.root", "RECREATE");

  TH1* sigClone = (TH1*) sig_->Clone(sig_->GetName()+TString("Stat"));
  DoStatVariation(sigClone, sig_->GetName());
  for(auto& bkg : bkg_)
    {
      TH1* bkgClone = (TH1*) bkg->Clone(bkg->GetName()+TString("Stat"));
      DoStatVariation(bkgClone, bkg->GetName());
    }
  data_->Write();
  sig_->Write();
  for(auto& bkg : bkg_)
    bkg->Write();
  for(auto& variation : statVariations_)
    variation->Write();
  
  shapesFile_->Close();
  if(debug_) cout << "[MakeSimpleCard::WriteShapes] Method finished successfully." << endl;
}

void MakeSimpleCard::DoStatVariation(TH1* shape, TString basename)
{
  ostringstream convert;   // stream used for bin name conversion
  vector<TH1*> variations;
  for(int ibin=1; ibin<=shape->GetNbinsX(); ++ibin)
    {
      convert << ibin;      // insert the textual representation of 'Number' in the characters in the stream
      TString binName(convert.str());
      TH1* tempUp   = (TH1*) shape->Clone(basename+TString("_bin")+binName+shape->GetName()+TString("Up"  ));
      TH1* tempDown = (TH1*) shape->Clone(basename+TString("_bin")+binName+shape->GetName()+TString("Down"));
      tempUp  ->SetBinContent(ibin, shape->GetBinContent(ibin)+shape->GetBinErrorUp(ibin));
      tempDown->SetBinContent(ibin, shape->GetBinContent(ibin)-shape->GetBinErrorLow(ibin));
      statVariations_.push_back(tempUp  );
      statVariations_.push_back(tempDown);
    }

}

void MakeSimpleCard::WriteCard()
{
  card_.close();
  cout << "[MakeSimpleCard::writeCard] Card has been written." << endl;
  cout << "                            Process it in CMSSW_7_1_X with:" << endl;
  cout << "                            combine -M Asymptotic " << cardName_ << ".txt" << endl;
  if(debug_) cout << "[MakeSimpleCard::WriteCard] Method finished successfully." << endl;
}
