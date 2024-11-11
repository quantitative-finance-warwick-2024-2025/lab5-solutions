#define CATCH_CONFIG_MAIN
#include "binomial_expectation.h"
#include "catch.hpp"
#include <vector>

// error tolerance for floating point comparisons
#define EPSILON 1e-14

void error(std::vector<double> &fxn, double p, unsigned int n) {
    double recursion_answer = binomial_expectation_recursion(fxn.data(), p, n);
    double induction_answer = binomial_expectation_loop(fxn.data(), p, n);

    // check answers are equal within the tolerance epsilon
    REQUIRE(induction_answer == Approx(recursion_answer).epsilon(EPSILON));
}

void populate_test_arrays(std::vector<double> &linear,
                          std::vector<double> &quadratic,
                          std::vector<double> &quartic,
                          unsigned int n) {
    for (unsigned int i = 0; i < n + 1; ++i) {
        linear[i] = (double)i - (double)n / 2.0;
        quadratic[i] = linear[i] * linear[i];
        quartic[i] = quadratic[i] * quadratic[i];
    }
}

TEST_CASE("recursion vs loop depth 1", "[tests]") {
    const unsigned int n = 1;
    double p = 0.5;
    std::vector<double> linear(n + 1);
    std::vector<double> quadratic(n + 1);
    std::vector<double> quartic(n + 1);

    populate_test_arrays(linear, quadratic, quartic, n);

    error(linear, p, n);
    error(quadratic, p, n);
    error(quartic, p, n);
}

TEST_CASE("recursion vs loop depth 3", "[tests]") {
    const unsigned int n = 3;
    double p = 0.5;
    std::vector<double> linear(n + 1);
    std::vector<double> quadratic(n + 1);
    std::vector<double> quartic(n + 1);

    populate_test_arrays(linear, quadratic, quartic, n);

    error(linear, p, n);
    error(quadratic, p, n);
    error(quartic, p, n);
}

TEST_CASE("Expectation and variance match theoretical values", "[analytic]") {
    const unsigned int n = 100;
    const double p = 0.9;

    // Calculate theoretical mean and variance
    double mean = n * p;
    double variance = n * p * (1.0 - p);

    // Create quadratic function (x - mean)^2
    std::vector<double> fx_x(n + 1);
    std::vector<double> fx_x_squared(n + 1);
    for (unsigned int i = 0; i < n + 1; ++i) {
        double x = static_cast<double>(i);
        fx_x[i] = x;
        fx_x_squared[i] = (x - mean) * (x - mean);
    }

    REQUIRE(binomial_expectation_loop(fx_x.data(), p, n) == Approx(mean).epsilon(EPSILON));
    // E[(X - mu)²] = E[X²] - (np)² =  [np(1-p) + (np)²] - (np)² = np(1-p)
    REQUIRE(binomial_expectation_loop(fx_x_squared.data(), p, n)
            == Approx(variance).epsilon(EPSILON));
}

void check_continuous_returns_error(double initial_value,
                                    unsigned int number_of_periods,
                                    double volatility,
                                    double time_period,
                                    double continuous_returns_rate) {
    double exact_value = std::pow(continuous_returns_rate, number_of_periods) * initial_value;
    double from_binomial_expectation = expected_returns(
        initial_value, number_of_periods, volatility, time_period, continuous_returns_rate);
    REQUIRE(from_binomial_expectation == Approx(exact_value).epsilon(EPSILON));
}

TEST_CASE("Continuous returns risk neutral rate volatility 1e0", "[analytic]") {
    double time_period = 4.0;
    unsigned int number_of_periods = 10;
    double risk_free_rate = 0.01;
    double continuous_returns_rate = std::exp(risk_free_rate * time_period / (double) number_of_periods);
    double initial_value = 100.0;
    check_continuous_returns_error(
        initial_value, number_of_periods, 1e0, time_period, continuous_returns_rate);
}

TEST_CASE("Continuous returns risk neutral rate volatility 1e1", "[analytic]") {
    double time_period = 4.0;
    unsigned int number_of_periods = 10;
    double risk_free_rate = 0.01;
    double continuous_returns_rate = std::exp(risk_free_rate * time_period / (double) number_of_periods);
    double initial_value = 100.0;
    check_continuous_returns_error(
        initial_value, number_of_periods, 1e1, time_period, continuous_returns_rate);
}

TEST_CASE("Continuous returns risk neutral rate volatility 1e2", "[analytic]") {
    double time_period = 4.0;
    unsigned int number_of_periods = 10;
    double risk_free_rate = 0.01;
    double continuous_returns_rate = std::exp(risk_free_rate * time_period / (double) number_of_periods);
    double initial_value = 100.0;
    check_continuous_returns_error(
        initial_value, number_of_periods, 1e2, time_period, continuous_returns_rate);
}