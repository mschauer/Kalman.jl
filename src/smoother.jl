
function smoother_kernel(Gs::T, Gf, Gpred, Φ) where {T}
    xs, Ps = meancov(Gs)
    xf, Pf = meancov(Gf)
    xpred, Ppred = meancov(Gpred)

    J = Pf*Φ'/Ppred # C/(C+w)
    xs = xf + J*(xs - xpred)
    Ps = Pf + J*(Ps - Ppred)*J'

    T(xs, Ps)
end

struct RauchTrungStriebel{T,S}
    iter::T # iterate backwards through
    M::S
end

function dyniterate(RTS::RauchTrungStriebel, ::Nothing)
    ϕ = dyniterate(RTS.iter, :Nothing)
    ϕ === nothing && return nothing
    (n, U), state = ϕ
    Gs, Gpred = U[1], U[2]
    n => Gs, (Gs, Gpred, state)
end


function dyniterate(RTS::RauchTrungStriebel, (Gs, Gpred, state))
    ϕ = dyniterate(RTS.iter, state)
    ϕ === nothing && return nothing
    (n, U), state = ϕ
    Gf, Gpredᵒ = U[1], U[2]
    Gs = smoother_kernel(Gs, Gf, Gpred, Phi(RTS.M.sys, n=>U))
    n => Gs, (Gs, Gpredᵒ, state)
end

function rts_smoother(M, prior, Y)
    P = filter(Y, M)

    ϕ = dyniterate(P, nothing, (value=prior,))
    ϕ === nothing && error("no observations")
    (t, u), state = ϕ

    Xf = trajectory((t => (u[1], u[2]),))
    while true
        ϕ = dyniterate(P, state)
        ϕ === nothing && break
        (t, u), state = ϕ
        push!(Xf, t => (u[1], u[2])) # safe only what is needed
    end
    ll = u[3]

    n = length(Xf)
    U = Xf.x[end]
    Gf, Gpred = U[1], U[2]
    Gs = Vector{typeof(Gf)}(undef, n)
    Gs[end] = Gf
    for i in n-1:-1:1
        U = Xf.x[i]
        Gf, Gpredᵒ = U[1], U[2]
        Gs[i] = smoother_kernel(Gs[i+1], Gf, Gpred, Phi(M.sys, i => U))
        Gpred = Gpredᵒ
    end
    Trajectory(Xf.t, Gs), ll
end
