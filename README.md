# IB9JHO Lab 5 - Binomial Expectation

## Key formulas:
- Time step: Δt = T/n
- Up factor: u = exp(σ√Δt)
- Down factor: d = 1/u
- Risk-neutral probability: p = (exp((r-q)Δt) - d)/(u - d)
- Terminal values: S₀ * u^(2j-n) for j ∈ [0,n]

## Exercises
Complete the following functions in `src/binomial_expectation.cpp`:

### 1. `binomial_expectation_loop`
Implement an iterative bottom-up computation of the binomial expectation.

**Recurrence relation:** V(i,j) = (1-p)V(i+1,j) + pV(i+1,j+1)  
**Terminal values:** V(n,j) = spread[j]

The implementation should work backwards from the terminal nodes to compute the root node value using an array to store intermediate values.

### 2. `binomial_expectation_recursion_h`
Implement a recursive top-down computation using the same recurrence relation and base case from step 1.
The implementation should recursively compute node values by traversing down to terminal nodes and combining results using the recurrence relation.

### 3. `expected_returns`
Implement the overloaded function that calculates values at spread_points number of points. For each point i, calculate:

- `out_spread[i]`: The spread factor from d to u, calculated as d + (u-d)*i/(spread_points-1)
- `out_risk_neutral_rate[i]`: The continuous rate corresponding to the spread, ln(spread[i])/Δt
- `out_probability[i]`: The probability corresponding to the spread, (spread[i]-d)/(u-d)
- `out_returns[i]`: The expected return using the spread as the continuous returns rate

## Testing and Submission
1. Run CTest to verify all tests pass
2. Push your changes with a commit message "Submit assignment"