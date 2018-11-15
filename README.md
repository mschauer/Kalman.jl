[![Build Status](https://travis-ci.org/mschauer/Kalman.jl.svg?branch=master)](https://travis-ci.org/mschauer/Kalman.jl)

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/Kalman.jl/latest/)


# Kalman
Flexible filtering and smoothing in Julia. `Kalman` uses [`DynamicIterators`](https://github.com/mschauer/DynamicIterators.jl) (an iterator protocol for dynamic data dependent and controlled processes) and
[`GaussianDistributions`](https://github.com/mschauer/GaussianDistributions.jl) (Gaussian distributions as abstraction for the uncertain state)
to implement flexible online Kalman filtering.

## Example

For the state space system

    x[k] = Φx[k−1] + b + w[k],    w[k] ∼ N(0, Q)

    y[k] = Hx[k] + v[k],    v[k] ∼ N(0, R)

define

```julia
using GaussianDistributions
using DynamicIterators

# Define linear evolution
Φ = [0.8 0.5; -0.1 0.8]
b = zeros(2)
Q = [0.2 0.0; 0.0 1.0]

E = LinearEvolution(Φ, Gaussian(b, Q))

# Define observation scheme
H = [1.0 0.0]
R = Matrix(1.0I, 1, 1)

O = LinearObservation(E, H, R)

# Prior
x0 = [1., 0.]
P0 = Matrix(1.0I, 2, 2)

# Observations (mock)
Y = [1 => [1.14326], 2 => [-0.271804], 3 => [-0.00512675]]

# Filter
Xf, ll = kalmanfilter(O, 0 => Gaussian(x0, P0), Y)
@show Xf

```

## Implementation
As said, filtering is implemented via the DynamicIterator protocol. It is worthwhile to look at
a possible the implementation of `kalmanfilter` to see how filtering can be integrated into online algorithms (run in a local scope to avoid `UndefVarError: ystate not defined`.)

```julia
# `Y` is the data iterator, iterating over pairs of  `t => v` of time `t` and observation `v`
# `O` is the dynamical filter iterator, iterating over pairs `t => u` where
#     u::Tuple{<:Gaussian,<:Gaussian,Float64}
# is the tuple of filtered state, the predicted state and the log likelihood

# Initialise data iterator

ϕ = iterate(Y)
ϕ === nothing && error("no observations")
(t, v), ystate = ϕ

# Initialise dynamical filter with first data point `t => v`
# and the `prior::Pair{Int,<:Gaussian}`, a pair of initial time and initial state

ϕ = dyniterate(O, Start(Kalman.Filter(prior, 0.0)), t => v)
ϕ === nothing && error("no observations")
(t, u), state = ϕ

X = trajectory((t => u[1],))
while true

    # Advance data iterator
    
    ϕ = iterate(Y, ystate)
    ϕ === nothing && break
    (t, v), ystate = ϕ

    # Advance filter with new data `t => v`
    
    ϕ = dyniterate(O, state, t => v)
    ϕ === nothing && break
    (t, u), state = ϕ
    
    # Do something with the result `t => u` (here: saving it)
    
    push!(X, t => u[1]) # save filtered state
end
ll = u[3] # likelihood
@show  X, ll
```
