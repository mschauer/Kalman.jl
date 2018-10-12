
struct GenericLinearObservation <: AbstractObservation
end
struct GenericLinearEvolution{TPhi,Tb,TQ} <: AbstractEvolution
end

"""
    LinearEvolution(Phi, b, Q)

Evolution `x -> Phi x + b + w` where ``w ~ N(0, Q)``
"""
struct LinearEvolution{TPhi,Tb,TQ} <: AbstractEvolution
    Phi::TPhi # dxd
    b::Tb # d
    Q::TQ # dxd
end

"""
    LinearObservation(H, R)

Observe `y = Hx + v` where ``v ~ N(0, R)``
"""
struct LinearObservation{TH, TR} <: AbstractObservation
    H::TH # dxd
    R::TR # d
end


function predict!(s, x, P, t, u, M) 
    if s == t # if `s == t`` no time has passed
        x + u, P, one(M.Phi)
    else
        M.Phi*x + M.b + u, M.Phi*P*M.Phi' + M.Q, M.Phi
    end
end

function predict!(s, x, P, t, U, M::GenericLinearEvolution) 
    Phi, b, Q = U # unpack generic input
    Phi*x + b, Phi*P*Phi' + Q, Phi
end

function observe!(s, x, P, t, Y, M::GenericLinearObservation)
    y, H, R = Y  # unpack generic observation
    t, y, H, R
end

function observe!(s, x, P, t, Y, MO::LinearObservation)
    t, Y, MO.H, MO.R
end

"""
    kalman_kernel(s, x, P, t, U, Y, SSM) -> t, x, P, Ppred, ll, K

Single Kalman filter step consisting of a prediction step `predict!`, an observation step `observe!`
and a correction step `correct!` with input `U`. Return filtered covariance `P` and predicted `Ppred`

Computes and returns as well the log likelihood of the residual and the Kalman gain.
"""
function kalman_kernel(s, x, P, t, U, Y, SSM)
   
    xpred, Ppred, Phi, b = predict!(s, x, P, t, U, SSM)

    t, y, H, R = observe!(s, xpred, P, t, Y, SSM)

    x, P, yres, S, K = correct!(xpred, Ppred, y, H, R, SSM)

    ll = llikelihood(yres, S, SSM)
    
    t, x, P, xpred, Ppred, ll, K
end

"""
    kalmanfilter!(tt, uu, yy, xx, PP, xxpred, PPpred, M)

Filter obserations in place.
`tt` are `n` (sic!) observation time points, `uu` are `n` control inputs and 
`yy` are `n` (generalized) observations (depending on context either vectors or matrices with `size(yy, 2) == n` etc.)
"""
function kalmanfilter!(tt, uu, yy, xx, PP, xxpred, PPpred, M) 
    n = length(tt)

    x = mean(prior(M))
    P = cov(prior(M))
    
    ll = 0.0  
    t = tt[1]

    for i in 1:n 
        t, x, P, xpred, Ppred, l, _ = kalman_kernel(t, x, P, tt[i], uu[.., i], yy[.., i], M)
        xx[.., i], PP[.., i], xxpred[.., i], PPpred[.., i] = x, P, xpred, Ppred                
        ll += l
    end
    xx, PP, xxpred, PPpred, ll
end


