
"""
    LinearEvolution(Φ, b, Q)

Evolution of the law of `x -> Φ x + w` where ``w ~ N(0, Q)``
"""
struct LinearEvolution{TΦ,TQ} <: Evolution
    Φ::TΦ # dxd
    Q::TQ # dxd
end

evolve(M::LinearEvolution, u::Pair) = timelift_evolve(M, u)

function evolve(M::LinearEvolution, G::Gaussian)
    Gaussian(M.Φ*mean(G) + mean(M.Q), M.Φ*cov(G)*M.Φ' + cov(M.Q))
end
