
struct End{T} <: Message
    state::T
end

struct Parent{T, S} <: Message
    state::T
    parent::S
end

struct Evolve{T} <: Message
    state::T
end

struct Trace{T,S,F} <: Message
    start::T
    X::S
    stop::F
end

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
