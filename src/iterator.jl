import Base: start, next, done, iteratorsize, iteratoreltype, length

struct KalmanFilter{Tit,TM}
    it::Tit
    M::TM
end

Base.iteratorsize(::KalmanFilter{Tit}) where {Tit} = Base.iteratorsize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit) 
Base.iteratoreltype(::KalmanFilter{Tit, LinearHomogSystem{Tx, TP}}) where {Tit, Tx, TP} = Gaussian{Tx,TP}

start(kf::KalmanFilter) = start(kf.it), prior(kf.M)
function next(kf::KalmanFilter, state)
    st, N = state

    y, st2 = next(kf.it, st)
    _, x, P = kalman_kernel(0.0, N.μ, N.σ, 1.0, y, SSM)
    Gaussian(x, P)
    (st2, N2), N2
end

done(kf::KalmanFilter, state) = done(kf.it, state[1])

length(kf::KalmanFilter) = length(kf)

struct TimedKalmanFilter{Tit,TM}
    it::Tit
    M::TM
end

Base.iteratorsize(::TimedKalmanFilter{Tit}) where {Tit} = Base.iteratorsize(Tit) == Base.HasShape() ? Base.HasLength() : Base.iteratorsize(Tit) 
Base.iteratoreltype(::TimedKalmanFilter{Tit, LinearHomogSystem{Tx, TP}}) where {Tit, Tx, TP} = Gaussian{Tx,TP}

start(kf::TimedKalmanFilter) = start(kf.it), -Inf, prior(kf.M)
function next(kf::TimedKalmanFilter, state)
    st, s, N = state

    ty, st2 = next(kf.it, st)
    t, y = ty
    t, x, P = kalman_kernel(s, N.μ, N.σ, t, y, SSM)
    Gaussian(x, P)
    (st2, N2), t, N2
end

done(kf::TimedKalmanFilter, state) = done(kf.it, state[1])

length(kf::TimedKalmanFilter) = length(kf)