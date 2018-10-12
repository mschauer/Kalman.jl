module Kalman

using Distributions, GaussianDistributions
export LinearHomogSystem, LinearStateSpaceModel, AbstractObservation, AbstractEvolution

export  kalmanfilter, kalmanfilter!, kalmanrts, kalmanrts!, kalmanEM

# generate
export sample, randmvn

# iterator
export KalmanFilter, MappedKalmanFilter

# track
export track

# input
export GenericLinearObservation, GenericLinearEvolution, LinearObservation, LinearEvolution

include("ellipse.jl")

abstract type FilterMethod
end

abstract type StateSpaceModel
end

abstract type AbstractObservation
end

abstract type AbstractEvolution
end

struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
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
struct LinearStateSpaceModel{Tsys,Tobs, Tprior} <: StateSpaceModel
    sys::Tsys
    obs::Tobs
    prior::Tprior
end
predict!(s, x, P, t, U, M::LinearStateSpaceModel) = predict!(s, x, P, t, U, M.sys)

predict!(s, x, P, t, M::LinearStateSpaceModel) = predict!(s, x, P, t, M.sys)
observe!(s, x, P, t, Y, M::LinearStateSpaceModel) = observe!(s, x, P, t, Y, M.obs)
prior(M::LinearStateSpaceModel) = M.prior

prior(M) = Gaussian(M.x0, M.P0)
dims(SSM) = size(SSM.H)
llikelihood(yres, S, SSM) = logpdf(Gaussian(zero(yres), S), yres)


function correct!(method::JosephForm, SSM, x, Ppred, y, H, R)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance
 
    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    x, P, yres, S, K
end


"""
    kalman_kernel(s, x, P, t, Y, SSM) -> t, x, P, Ppred, ll, K

Single Kalman filter step consisting of a prediction step `predict!`, an observation step `observe!`
and a correction step `correct!`. Return filtered covariance `P` and predicted `Ppred`

Computes and returns as well the log likelihood of the residual and the Kalman gain.
"""
function kalman_kernel(s, x, P, t, Y, SSM)
   
    x, Ppred, Phi = predict!(s, x, P, t, SSM)

    t, y, H, R = observe!(s, x, P, t, Y, SSM)

    x, P, yres, S, K = correct!(x, Ppred, y, H, R, SSM)

    ll = llikelihood(yres, S, SSM)
    
    t, x, P, Ppred, ll, K
end


function smoother_kernel(xs, Ps, xf, Pf, Ppred, Phi, b)
    # xs -- previous #h[i+1]
    # Ps -- previous #H[i+1]
    # xf -- xxf[:, i] #m[i]
    # Pf -- PPf[:, :, i] # Cn
    # Ppred -- PPpred[:, :, i+1] # R[i+1] = C[i] + W[i+1]

    J = Pf*Phi'/Ppred # C/(C+w)
    xs = xf +  J*(xs - (Phi*xf  + b))
    Ps = Pf + J*(Ps - Ppred)*J' 
   
    xs, Ps, J
end

function backward_kernel(x, Ps, xf, Pf, Ppred, Phi, b)
    # x -- previous
    # Ps -- previous
    # xf -- xxf[:, i]
    # Pf -- PPf[:, :, i]
    # Ppred -- PPpred[:, :, i+1]

    J = Pf*Phi'/Ppred
    xs = xf +  J*(x - (Phi*xf  + b))
    Ps = Pf - J*Ppred*J' # as Ppred = M.Phi*P*M.Phi' + M.Q
    xs, Ps, J
end




include("input.jl")

include("linearhomogsystem.jl")

include("generate.jl")

include("iterator.jl")

include("track.jl")

#include("kalmanem.jl")


end # module
