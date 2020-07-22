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

#include "CustomPredictionStrategy.h"

namespace grf {

size_t CustomPredictionStrategy::prediction_length() const {
  return 1;
}

std::vector<double> CustomPredictionStrategy::predict(size_t sample,
    const std::unordered_map<size_t, double>& weights_by_sample,
    const Data& train_data,
    const Data& data) const {

  size_t num_samples = weights_by_sample.size();
  double total_treatment = 0.0;
  size_t number_of_outcomes = 2; // generalize // you CAN have custom prediction length // or query
  double total_effect = 0.0;
  std::vector<double> effects(number_of_outcomes);

  for (const auto& entry : weights_by_sample) {
    total_treatment += data.get_treatment(sample);
  }

  for (const auto& entry : weights_by_sample) {
    size_t sample = entry.first;
    double denominator = (
      data.get_treatment(sample) * total_treatment -
      (num_samples - total_treatment) * (1 - data.get_treatment(sample))
    );
    effects[(size_t) data.get_instrument(sample)] += data.get_outcome(sample) / denominator;
    total_effect += data.get_outcome(sample) / denominator;
  }
  return { effects[0] / total_effect };
}

std::vector<double> CustomPredictionStrategy::compute_variance(
    size_t sample,
    const std::vector<std::vector<size_t>>& samples_by_tree,
    const std::unordered_map<size_t, double>& weights_by_sampleID,
    const Data& train_data,
    const Data& data,
    size_t ci_group_size) const {
  return { 0.0 };
}

} // namespace grf
