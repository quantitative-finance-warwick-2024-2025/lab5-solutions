#include "binomial_expectation.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>

void write_binomial_results(const std::filesystem::path &project_dir,
                            const std::string &filename,
                            double S0,
                            unsigned int number_of_periods,
                            double volatility,
                            double time_period,
                            size_t spread_points) {
    // Create output directory
    std::filesystem::path output_dir = project_dir / "output";
    std::filesystem::create_directories(output_dir);

    // Prepare output file - add .csv extension if not present
    std::string final_filename = filename;
    if (!filename.ends_with(".csv")) {
        final_filename += ".csv";
    }
    std::filesystem::path output_file = output_dir / final_filename;
    std::ofstream file(output_file);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open output file: " + output_file.string());
    }

    // Allocate arrays for results
    std::vector<double> out_spread(spread_points);
    std::vector<double> out_returns(spread_points);
    std::vector<double> out_probability(spread_points);
    std::vector<double> out_risk_neutral_rate(spread_points);

    // Calculate expected returns for different spread points
    expected_returns(S0,
                     number_of_periods,
                     volatility,
                     time_period,
                     spread_points,
                     out_spread.data(),
                     out_risk_neutral_rate.data(),
                     out_probability.data(),
                     out_returns.data());

    // Write CSV header
    file << "out_spread,risk_neutral_rate,probability,log_returns_over_period\n";

    // Write results with high precision
    file << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < spread_points; ++i) {
        file << out_spread[i] << "," << out_risk_neutral_rate[i] << "," << out_probability[i] << ","
             << std::log(out_returns[i]) << "\n";
    }

    file.close();
}

int main() {
    std::filesystem::path project_dir(PROJECT_DIR);

    // Example with different filenames
    write_binomial_results(project_dir,
                           "binomial_results_10periods", // filename without extension
                           100.0,                         // S0
                           10,                           // number_of_periods
                           0.5,                           // volatility
                           5.0,                          // time_period
                           10);                          // spread_points

    write_binomial_results(project_dir,
                           "binomial_results_100periods.csv", // filename with extension
                           100.0,                            // S0
                           100,                               // number_of_periods
                           0.5,                              // volatility
                           5.0,                              // time_period
                           10);                             // spread_points

    write_binomial_results(project_dir,
                           "binomial_results_1000periods.csv", // filename with extension
                           100.0,                            // S0
                           1000,                               // number_of_periods
                           0.5,                              // volatility
                           5.0,                              // time_period
                           10);                             // spread_points

    return 0;
}