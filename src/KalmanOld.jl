module Kalman

using Distributions, GaussianDistributions
using Trajectories, DynamicIterators
using Random, LinearAlgebra

import Distributions: sample

export LinearHomogSystem, LinearStateSpaceModel, AbstractObservation, AbstractEvolution

export kalmanfilter, kalmanfilter!, kalmanrts, kalmanrts!, kalmanEM

# generate
export sample, randmvn

# iterator
export KalmanFilter, MappedKalmanFilter

# track
export track

# input
export GenericLinearObservation2, GenericLinearEvolution, LinearObservation2, LinearEvolution

include("ellipse.jl")

meancov(G) = mean(G), cov(G)
meancov(G::Tuple) = G



abstract type StateSpaceModel
end


abstract type AbstractEvolution
end
abstract type FilterMethod
end
struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
end
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





include("general.jl")


include("iterator.jl")

include("track.jl")

#include("kalmanem.jl")

end # module
