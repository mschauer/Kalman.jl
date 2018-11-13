
struct End{T} <: Message
    state::T
end

struct Parent{T, S} <: Message
    state::T
    parent::S
end
struct Sample{T,RNG} <: Message
    state::T
    rng::RNG
end

Sample(x) = Sample(x, Random.GLOBAL_RNG)

struct Evolve{T} <: Message
    state::T
end

struct Trace{T,S,F} <: Message
    start::T
    X::S
    stop::F
end

iterate(M::Message) = getfield(M, 1), 1
iterate(M::Message, n) = getfield(M, n + 1), n + 1

import DynamicIterators: trace
function trace(iter, (start, X, stop)::Trace)
    ϕ = dyniterate(iter, start)
    ϕ === nothing && return X
    x, u = ϕ
    push!(X, x)
    while !stop(x)
        ϕ = dyniterate(iter, u)
        ϕ === nothing && break
        x, u = ϕ
        push!(X, x)
    end
    return X
end

export Start, Trace, Sample
