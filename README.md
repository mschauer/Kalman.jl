[![Build Status](https://travis-ci.org/mschauer/Kalman.jl.svg?branch=master)](https://travis-ci.org/mschauer/Kalman.jl)

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/Kalman.jl/latest/)


# Kalman

`Kalman.jl` uses `DynamicIterators` to implement online Kalman filtering.

## Example

For the state space system

    x[k] = Φx[k−1] + b + w[k],    w[k] ∼ N(0, Q)

    y[k] = Hx[k] + v[k],    v[k] ∼ N(0, R)

define

```julia
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
ϕ = iterate(Y)
ϕ === nothing && error("no observations")
y, ystate = ϕ

ϕ = dyniterate(O, Start(Kalman.Filter(prior, 0.0)), y)
ϕ === nothing && error("no observations")
(t, u), state = ϕ

X = trajectory((t => u[1],))
while true
    ϕ = iterate(Y, ystate)
    ϕ === nothing && break
    y, ystate = ϕ

    ϕ = dyniterate(O, state, y)
    ϕ === nothing && break
    (t, u), state = ϕ
    push!(X, t => u[1]) # filtered state as Gaussian
                        # the second argument is predicted state
end
ll = u[3] # likelihood
@show  X, ll

```
