
abstract type FilterMethod
end
struct JosephForm <: FilterMethod
end
struct SimpleKalman <: FilterMethod
end

function correct(SSM::LinearStateSpaceModel, u::T, y) where T
    x, Ppred = meancov(u)
    H, R = SSM.obs.H, SSM.obs.R
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S
end

correct(O::LinearObservation, u, y) = correct(JosephForm(), u, (Gaussian(y, O.R), O.H))

function correct(method::JosephForm, u::T, (v, H)::Tuple{<:Gaussian, <:Any}) where T
    x, Ppred = meancov(u)
    y, R = meancov(v)
    yres = y - H*x # innovation residual

    S = H*Ppred*H' + R # innovation covariance

    K = Ppred*H'/S # Kalman gain
    x = x + K*yres
    P = (I - K*H)*Ppred*(I - K*H)' + K*R*K'
    T(x, P), yres, S
end

"""

Single Kalman filter with control being the observation, consisting of a prediction step
(evolution of the underlying system) and a correction step `correct`.
Return filtered state `t => (x, P)` and predicted `(xpred, Ppred)`.

Computes and returns as well the log likelihood of the residual.
"""
function dyniterate(M::StateSpaceModel, ((s, u), ll)::Condition{<:Pair}, v)
    t, y = v
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S = correct(M, upred, y)
    llᵒ = llikelihood(M, yres, S)
    (t => (u, upred, ll + llᵒ)), Condition(t => u, ll + llᵒ)
end

function dyniterate(M::StateSpaceModel, ((value,), dummy)::Condition{<:Start}, v)
    s, u = value
    t, y = v
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S = correct(M, upred, y)
    ll = llikelihood(M, yres, S)
    (t => (u, upred, ll)), Condition(t => u, ll)
end



function dyniterate(O::LinearObservation{<:LinearEvolution}, ((s, u), ll)::Condition, v)
    t, y = v
    t, upred = evolve(O.P, s => u, t)
    u, yres, S = correct(O, upred, y)
    llᵒ = llikelihood(O, yres, S)
    (t => u), Condition(t => u, ll + llᵒ)
end


function dyniterate(O::LinearObservation{<:LinearEvolution}, ((s, u), ll)::Filter, v)
    t, y = v
    t, upred = evolve(O.P, s => u, t)
    u, yres, S = correct(O, upred, y)
    llᵒ = llikelihood(O, yres, S)
    (t => (u, upred, ll + llᵒ)), Filter(t => u, ll + llᵒ)
end

function dyniterate(O::LinearObservation{<:LinearEvolution}, start::Start{<:Filter}, v)
    ((s, u), ll) = start.value
    t, y = v
    t, upred = evolve(O.P, s => u, t)
    u, yres, S = correct(O, upred, y)
    llᵒ = llikelihood(O, yres, S)
    (t => (u, upred, ll + llᵒ)), Filter(t => u, ll + llᵒ)
end

llikelihood(yres, S) = logpdf(Gaussian(zero(yres), symmetrize!(S)), yres)
llikelihood(_, yres, S) = logpdf(Gaussian(zero(yres), symmetrize!(S)), yres)


#=
function dyniterate(M::StateSpaceModel, (s, u)::Pair, (v,)::Observe)
    t, (v, H) = v
    t, upred = evolve(M.sys, s => u, t)
    u, yres, S, K = correct(JosephForm(), upred, (v, H))
    ll = llikelihood(M, yres, S)
    t => (u, Ppred, ll, K), t => u
end
=#


"""
    kalmanfilter(M, t => x0) -> kf

`kf(iter)` is an iterator over `Gaussian`s or `Distributions` representing
the filtered distribution of `x` where `y` iterates over
(enumerated) signal values.

# Example

```
kf = kalmanfilter(M, 0 => prior) #
est1 = collect(kf(Y1))
est2 = collect(kf(Y2))
```
"""
kalmanfilter(M, prior) = iter -> filter(iter, from(M, prior))
kalmanfilter(M) = iter -> filter(iter, M)

function kalmanfilter(M, prior, Y)
    P = filter(Y, M)

    ϕ = dyniterate(P, Start(prior))
    ϕ === nothing && error("no observations")
    (t, u), state = ϕ

    X = trajectory((t => u[1],))
    while true
        ϕ = dyniterate(P, state)
        ϕ === nothing && break
        (t, u), state = ϕ
        push!(X, t => u[1])
    end
    ll = u[3]
    X, ll
end

function kalmanfilter(O::Observation, prior, Y)
    P = Bind(Y, O)

    ϕ = dyniterate(P, Start(Filter(prior, 0.0)))
    ϕ === nothing && error("no observations")
    (t, u), state = ϕ

    X = trajectory((t => u[1],))
    while true
        ϕ = dyniterate(P, state)
        ϕ === nothing && break
        (t, u), state = ϕ
        push!(X, t => u[1])
    end
    ll = u[3]
    X, ll
end
