#include <cmath>
#include <vector>

double binomial_expectation_loop(double *spread, double p, unsigned int n) {
    std::vector<double> arr(n + 1);

    // set up the base nodes
    for (unsigned int i = 0; i < n + 1; ++i)
        arr[i] = spread[i];

    // loop through each depth
    for (int i = (int)n - 1; i >= 0; --i)
        // loop through each node at depth i
        for (unsigned int j = 0; j < (unsigned int)i + 1; ++j)
            arr[j] = (1.0 - p) * arr[j] + p * arr[j + 1];

    double answer = arr[0];
    return answer;
}

double binomial_expectation_recursion_h(
    unsigned int i, unsigned int j, double *spread, double p, unsigned int n) {
    if (i == n)
        return spread[j];

    // if (i < n)
    return (1.0 - p) * binomial_expectation_recursion_h(i + 1, j, spread, p, n)
           + p * binomial_expectation_recursion_h(i + 1, j + 1, spread, p, n);
}

double binomial_expectation_recursion(double *spread, double p, unsigned int n) {
    return binomial_expectation_recursion_h(0, 0, spread, p, n);
}

double risk_neutral_probability(double volatility, double delta_t, double rate_of_return) {
    double u = std::exp(volatility * std::sqrt(delta_t));
    double d = 1.0 / u;
    return (std::exp(rate_of_return) - d) / (u - d);
}

double expected_returns(double S0,
                        unsigned int number_of_periods,
                        double volatility,
                        double time_period,
                        double continuous_returns_rate) {
    double delta_t = time_period / (double)number_of_periods;
    double u = std::exp(volatility * std::sqrt(delta_t));
    double d = 1.0 / u;
    double p = (continuous_returns_rate - d) / (u - d);

    std::vector<double> spread(number_of_periods + 1);

    for (unsigned int i = 0; i < number_of_periods + 1; ++i) {
        spread[i] = S0 * std::pow(u, 2.0 * (double)i - (double)number_of_periods);
    }

    return binomial_expectation_loop(spread.data(), p, number_of_periods);
}

void expected_returns(double S0,
                      unsigned int number_of_periods,
                      double volatility,
                      double time_period,
                      size_t spread_points,
                      double *out_spread,
                      double *out_risk_neutral_rate,
                      double *out_probability,
                      double *out_returns) {

    double delta_t = time_period / (double)number_of_periods;
    double u = std::exp(volatility * std::sqrt(delta_t));
    double d = 1.0 / u;

    for (size_t i = 0; i < spread_points; ++i) {
        out_spread[i] = d + (u - d) * i / (spread_points - 1);
        out_risk_neutral_rate[i] = std::log(out_spread[i]) / delta_t;
        out_probability[i] = (out_spread[i] - d) / (u - d);
        out_returns[i] =
            expected_returns(S0, number_of_periods, volatility, time_period, out_spread[i]);
    }
}
