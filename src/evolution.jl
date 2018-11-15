
"""
    LinearEvolution(Φ, b, Q) <: Evolution

Evolution of the law of `x -> Φ x + w` where ``w \\sim N(0, Q)``.

# Examples
```
    evolve(LinearEvolution(Φ, b, Q), 0 => Gaussian(x, P))
```
"""
struct LinearEvolution{TΦ,TQ} <: Evolution
    Φ::TΦ # dxd
    Q::TQ # dxd
end
Phi(M::LinearEvolution, ::Any) = M.Φ
Phi(M::LinearEvolution) = M.Φ

struct GenericLinearEvolution <: Evolution
end

evolve(M::LinearEvolution, u::Pair) = timelift_evolve(M, u)
evolve(M::LinearEvolution, u::Pair, c) = timelift_evolve(M, u, c)

function evolve(M::LinearEvolution, G)
    M.Φ*G ⊕ M.Q
end


function evolve(M::GenericLinearEvolution, (G,)::Control, control)
    Φ, Q = control
    G(Φ*mean(G) + mean(Q), Φ*cov(G)*Φ' + cov(Q))
end
