#include "Significance.h"
#include <iostream>
#include <TString.h>
#include <TMath.h>

Significance::Significance(std::vector<double> input):
  significance_(-9.),
  S_(-9.),
  B_(-9.),
  dB_(0.)
{
  SetYields(input);

}

void Significance::SetYields(std::vector<double> input)
{
  input_=input;
  if(input_.size()<2)
    cout << "[Significance::SetYields] ERROR: I don't have both signal and background. Input size " << input_.size() << endl;
  else if(input_.size()>=2)
    {
      S_=input_[0];
      B_=input_[1];
    }
  if(input_.size()>2)
    dB_=input_[2];
  
  cout << "[Significance::SetYields] Signal events: " << S_ << ", Background events: " << B_ << endl;  
  
}


double Significance::getSignificance(TString which, int a, int b)
{

  if(input_.size()<2)
    {
      cout << "[Significance::Significance] ERROR: I don't have both signal and background. Input size " << input_.size() << endl;
      return -9.;
    }
  if(which=="punzi")
    ComputePunzi(a, b);
  else if(which=="approxpunzi")
    ComputeApproximatedPunzi(a);
  else if(which=="pseudosearch")
    ComputePseudoSignificanceSearchOriented();
  else if(which=="pseudodiscovery")
    ComputePseudoSignificanceDiscoveryOriented();
  else if (which=="pvalue")
    ComputePValue();
  else
    cout << "Please select an expression for the significance" << endl;
  
  return significance_;
}


void Significance::ComputePunzi(int a, int b)
{
  // Return eq.6 of arxiv:physics/0308063v2
  significance_ = S_/( b*b + 2*a*sqrt(B_) + b*sqrt(b*b + 4*a*sqrt(B_) + 4*B_));
}

void Significance::ComputeApproximatedPunzi(int a)
{
  // Return eq.7 of arxiv:physics/0308063v2
  significance_ = S_/(a/2 + sqrt(B_));

}

void Significance::ComputePValue()
{
  significance_ = TMath::Poisson( S_ + B_, B_);
}

void Significance::ComputePseudoSignificanceSearchOriented()
{
  // Return s/sqrt(b)
  significance_ = S_/sqrt(B_);
}

void Significance::ComputePseudoSignificanceDiscoveryOriented()
{
  // Return s/sqrt(s+b)
  significance_ = S_/(sqrt(S_ + B_));
}

void Significance::Test()
{
  std::vector<double> i1 = {23.};
  std::vector<double> i2 = {24., 4560.};
  std::vector<double> i3 = {25., 4560., 123.};

  SetYields(i1);
  std::cout << "1d input array" << std::endl;
  std::cout << "\t punzi:           " << getSignificance("punzi")           << endl;
  std::cout << "\t approxpunzi:     " << getSignificance("approxpunzi")     << endl;
  std::cout << "\t pseudosearch:    " << getSignificance("pseudosearch")    << endl;
  std::cout << "\t pseudodiscovery: " << getSignificance("pseudodiscovery") << endl;

  SetYields(i2);
  std::cout << "2d input array" << std::endl;
  std::cout << "\t punzi:           " << getSignificance("punzi")           << endl;
  std::cout << "\t approxpunzi:     " << getSignificance("approxpunzi")     << endl;
  std::cout << "\t pseudosearch:    " << getSignificance("pseudosearch")    << endl;
  std::cout << "\t pseudodiscovery: " << getSignificance("pseudodiscovery") << endl;

  SetYields(i3);
  std::cout << "3d input array" << std::endl;
  std::cout << "\t punzi:           " << getSignificance("punzi")           << endl;
  std::cout << "\t approxpunzi:     " << getSignificance("approxpunzi")     << endl;
  std::cout << "\t pseudosearch:    " << getSignificance("pseudosearch")    << endl;
  std::cout << "\t pseudodiscovery: " << getSignificance("pseudodiscovery") << endl;

}
