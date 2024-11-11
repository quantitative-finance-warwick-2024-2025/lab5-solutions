#pragma once

/**
 * @brief Calculates expected value in a binomial tree using an iterative bottom-up approach
 * @param spread Array of terminal node values (size n+1)
 * @param p Probability of up movement
 * @param n Number of periods in binomial tree
 * @return Expected value at root node
 */
double binomial_expectation_loop(double *spread, double p, unsigned int n);

/**
 * @brief Calculates expected value in a binomial tree using recursion
 * @param spread Array of terminal node values (size n+1)
 * @param p Probability of up movement
 * @param n Number of periods in binomial tree
 * @return Expected value at root node
 */
double binomial_expectation_recursion(double *spread, double p, unsigned int n);

/**
 * @brief Calculates risk-neutral probability in a binomial model
 * @param volatility Asset price volatility
 * @param delta_t Time step size
 * @param rate_of_return Continuous rate of return
 * @return Risk-neutral probability of up movement
 */
double risk_neutral_probability(double volatility, double delta_t, double rate_of_return);

/**
 * @brief Calculates expected asset price using binomial model
 * @param S0 Initial asset price
 * @param number_of_periods Number of time periods
 * @param volatility Asset price volatility
 * @param time_period Total time period length
 * @param continuous_returns_rate Continuous rate of return
 * @return Expected asset price
 */
double expected_returns(double S0,
                       unsigned int number_of_periods,
                       double volatility,
                       double time_period,
                       double continuous_returns_rate);

/**
 * @brief Calculates spread of expected returns and associated metrics
 * @param S0 Initial asset price
 * @param number_of_periods Number of time periods
 * @param volatility Asset price volatility
 * @param time_period Total time period length
 * @param spread_points Number of points in spread calculation
 * @param out_spread Output array for spread factors from d to u
 * @param out_risk_neutral_rate Output array for corresponding continuous rates
 * @param out_probability Output array for corresponding probabilities
 * @param out_returns Output array for expected returns at each spread point
 */
void expected_returns(double S0,
                     unsigned int number_of_periods,
                     double volatility,
                     double time_period,
                     size_t spread_points,
                     double *out_spread,
                     double *out_risk_neutral_rate,
                     double *out_probability,
                     double *out_returns);