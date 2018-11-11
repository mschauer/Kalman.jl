abstract type FilterMethod
end
struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
end

function correct!(method::JosephForm, SSM, Gpred::T, y, H, R) where T
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
