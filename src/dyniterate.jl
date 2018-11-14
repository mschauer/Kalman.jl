

struct Trace{T,S,F} <: Message
    start::T
    X::S
    stop::F
end

function collecttrace(iter, t::Trace)
    start = t.start
    X = t.X
    stop = t.stop
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

export Trace, collecttrace
