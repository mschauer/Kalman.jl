[![Build Status](https://travis-ci.org/mschauer/Kalman.jl.svg?branch=master)](https://travis-ci.org/mschauer/Kalman.jl)

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/Kalman.jl/latest/)


# Kalman

## Example

For the state space system

    x[k] = Φx[k−1] + b + w[k],    w[k] ∼ N(0, Q)
 
    y[k] = Hx[k] + v[k],    v[k] ∼ N(0, R)

define

```
x0 = [1., 0.]
P0 = eye(2)

Phi = [0.8 0.2; 0.0 0.8]
b = zeros(2)
Q = [0.1 0.0; 0.0 1.0]

y = [NaN]
H = [1.0 0.0]
R = eye(1)
M = LinearHomogSystem(x0, P0, Phi, b, Q, y, H, R)
```

and filter vector of observations `Y = [y[k] for k in 1:n]`  as

```
kf = KalmanFilter(Y, M)

est = collect(kf)
```

## Installation

```
Pkg.clone("https://github.com/mschauer/GaussianDistributions.jl")
Pkg.clone("https://github.com/mschauer/Kalman.jl")
```