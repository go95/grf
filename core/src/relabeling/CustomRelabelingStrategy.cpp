/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/
#include <vector>
#include <cmath>

#include "CustomRelabelingStrategy.h"

namespace grf {

// zero division errors etc!
// weight

bool CustomRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    std::vector<double>& responses_by_sample) const {

  // Change instrument to proper cluster
  size_t num_samples = samples.size();
  size_t number_of_outcomes = 2; // generalize // you CAN have custom prediction length // or query
  double total_treatment = 0.0;
  double total_effect = 0.0;
  std::vector<double> effects(number_of_outcomes);

  for (size_t i = 0; i < number_of_outcomes; i++) {
    total_outcome[i] = 0;
  }

  for (size_t sample : samples) {
    total_treatment += data.get_treatment(sample);
  }

  for (size_t sample : samples) {
    double denominator = (
      data.get_treatment(sample) * total_treatment -
      (num_samples - total_treatment) * (1 - data.get_treatment(sample))
    );
    effects[(size_t) data.get_instrument(sample)] += data.get_outcome(sample) / denominator;
    total_effect += data.get_outcome(sample) / denominator;
  }

  for (size_t sample : samples) {
    double denominator = (
      data.get_treatment(sample) * total_treatment -
      (num_samples - total_treatment) * (1 - data.get_treatment(sample))
    );

    responses_by_sample[sample] = (
      (1 - effects[data.get_instrument(sample)] / total_effect) * 
      data.get_outcome(sample) / denominator
    );
  }
  return false;
}

} // namespace grf
