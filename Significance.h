#ifndef SIGNIFICANCE_H
#define SIGNIFICANCE_H

// Forward declarations
class TString;

class Significance
{
 public:
  Significance(std::vector<double> input={});
  ~Significance() {};
  
  double getSignificance(TString which, int a=2, int b=2);
  void Test();

 protected:
  
  std::vector<double> input_;
  double S_;
  double B_;
  double dB_;

  double significance_;

  void SetYields(std::vector<double> input);
  void ComputePunzi(int a, int b);
  void ComputeApproximatedPunzi(int a);
  void ComputePseudoSignificanceSearchOriented();
  void ComputePseudoSignificanceDiscoveryOriented();
};

#endif

