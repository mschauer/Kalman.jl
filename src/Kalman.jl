module Kalman

using Distributions, GaussianDistributions
using Trajectories, DynamicIterators
using Random, LinearAlgebra

import DynamicIterators: evolve, dyniterate
import Base: iterate, IteratorSize, IteratorEltype, eltype, length

import Random.rand

meancov(G) = mean(G), cov(G)
meancov(G::Tuple) = G

export LinearObservation, GenericLinearObservation, Observation
include("observation.jl")

export LinearEvolution, Evolution
include("evolution.jl")

export LinearStateSpaceModel, StateObs
include("statespacemodel.jl")

include("filter.jl")
include("smoother.jl")
include("backwardsampler.jl")


#=
import Distributions: sample

export LinearHomogSystem, LinearStateSpaceModel, AbstractObservation, AbstractEvolution

export kalmanfilter, kalmanfilter!, kalmanrts, kalmanrts!, kalmanEM

# generate
export sample, randmvn

# iterator
export KalmanFilter, MappedKalmanFilter

# input
export GenericLinearObservation, GenericLinearEvolution, LinearObservation, LinearEvolution

include("ellipse.jl")



abstract type StateSpaceModel
end

"""
```
LinearStateSpaceModel <: StateSpaceModel

LinearStateSpaceModel(sys, obs, prior)
```

Combines a linear system `sys`, an observations model `obs` and
a `prior` to a linear statespace model in a modular way. See [LinearHomogSystem`](@ref)
for a "batteries included" complete linear system.
"""
struct LinearStateSpaceModel{T,Tsys,Tobs} <: StateSpaceModel
    prior::T
    sys::Tsys
    obs::Tobs
end

prior(M) = M.prior
gausstype(::LinearStateSpaceModel{T}) where {T} = T

predict!(s, G, t, U, M::LinearStateSpaceModel) = predict!(s, G, t, U, M.sys)

predict!(s, G, t, M::LinearStateSpaceModel) = predict!(s, G, t, M.sys)
observe!(s, t, Y, M::LinearStateSpaceModel) = observe!(s, t, Y, M.obs)
prior(M::LinearStateSpaceModel) = M.prior

dims(SSM) = size(SSM.H)
llikelihood(yres, S, SSM) = logpdf(Gaussian(zero(yres), S), yres)


function correct!(method::JosephForm, SSM, Gpred::T, y, H, R) where T
    x, Ppred = meancov(Gpred)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S, K
end

function correct!(method::JosephForm, SSM::LinearStateSpaceModel{T}, Gpred, y, H, R) where T
    x, Ppred = meancov(Gpred)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S, K
end

"""
    kalman_kernel(s, G, t, Y, SSM) -> t, G, Ppred, ll, K

Single Kalman filter step consisting of a prediction step `predict!`, an observation step `observe!`
and a correction step `correct!`. Return filtered covariance `P` and predicted `Ppred`

Computes and returns as well the log likelihood of the residual and the Kalman gain.
"""
function kalman_kernel(s, G, t, Y, SSM)

    Gpred, Phi = predict!(s, G, t, SSM)

    t, y, H, R = observe!(s, t, Y, SSM)

    G, yres, S, K = correct!(Gpred, y, H, R, SSM)

    ll = llikelihood(yres, S, SSM)

    t, G, Ppred, ll, K
end


function smoother_kernel(Gs::T, Gf, Ppred, Phi, b) where {T}

    # xs -- previous #h[i+1]
    # Ps -- previous #H[i+1]
    # xf -- xxf[:, i] #m[i]
    # Pf -- PPf[:, :, i] # Cn
    # Ppred -- PPpred[:, :, i+1] # R[i+1] = C[i] + W[i+1]
    xs, Ps = meancov(Gs)
    xf, Pf = meancov(Gf)

    J = Pf*Phi'/Ppred # C/(C+w)
    xs = xf +  J*(xs - (Phi*xf  + b))
    Ps = Pf + J*(Ps - Ppred)*J'

    T(xs, Ps), J
end

function backward_kernel(Gs::T, Gf, Ppred, Phi, b) where T
    # x -- previous
    # Ps -- previous
    # xf -- xxf[:, i]
    # Pf -- PPf[:, :, i]
    # Ppred -- PPpred[:, :, i+1]
    xs, Ps = meancov(Gs)
    xf, Pf = meancov(Gf)

    J = Pf*Phi'/Ppred
    xs = xf +  J*(x - (Phi*xf  + b))
    Ps = Pf - J*Ppred*J' # as Ppred = M.Phi*P*M.Phi' + M.Q
    T(xs, Ps), J
end




include("general.jl")


include("iterator.jl")

include("track.jl")

#include("kalmanem.jl")
=#
end # module
