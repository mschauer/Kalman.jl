
function backward_rand_kernel(x, Gf, Gpred, Phi) where T
    xf, Pf = meancov(Gf)
    xpred, Ppred = meancov(Gpred)

    J = Pf*Phi'/Ppred
    xs = xf +  J*(x - (Phi*xf  + b))
    Ps = Pf - J*Ppred*J' # as Ppred = M.Phi*P*M.Phi' + M.Q
    rand(Gaussian(xs, Ps))
end

struct BackwardRand{T,S}
    iter::T # iterate backwards through
    M::S
end

function dyniterate(B::BackwardRand, (x, state))
    ϕ = dyniterate(B.iter, state)
    ϕ === nothing && return nothing
    (n, U), state = ϕ
    Gf, Gpred = U[1], U[2]
    Gs = backward_rand_kernel(x, Gf, Gpred, Phi(RTS.M.sys, (n, U)))
    n=>x, (x, state)
end
