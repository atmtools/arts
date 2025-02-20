
#include "atm.h"
#include "psd.h"

#include "test_utils.h"

bool test_binned_psd() {

  Vector bins = {0.1, 0.5, 1.0};
  Vector counts = {1.0, 2.0};

  scattering::BinnedPSD psd{SizeParameter::DVeq, bins, counts, 272.5, 350};

  auto point = AtmPoint{1e4, 300};

  Vector sizes = {0.05, 0.3, 0.75, 1.1};
  Vector pnd = psd.evaluate(point, sizes, 1.0, 1.0);
  Vector pnd_ref = {0.0, 1.0, 2.0, 0.0};

  auto err = max_error(pnd, pnd_ref);
  if (err > 1e-6) return false;
  return true;
}


int main() {
  bool passed = false;
  std::cout << "Test Binned PSD:\t";
  passed = test_binned_psd();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
  return 0;
}
