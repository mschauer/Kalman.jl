
"""
    filter(Y, P)

"Filter" data `Y` with iterator `P`
calling (handling of `nothing`s omitted)
```
(t, y), state = iterate(Y, state)
(s => x) = dyniterate(P, s => x, (observation = t => y,))
```

"""
struct Filtered{T,S} <: DynamicIterator
    Y::S
    P::T
end
filter(Y, P, ::Nothing) = Filtered(Y, P)

function dyniterate(M::Filtered, (value,)::Start)
    v, q = @returnnothing dyniterate(M.Y)
    u, p = @returnnothing dyniterate(M.P, nothing, (value = value, observation = v))
    u, (q, p)
end

function iterate(M::Filtered, (q, p)::Tuple)
    v, q = @returnnothing iterate(M.Y, q)
    u, p = @returnnothing dyniterate(M.P, p, (observation = v,))
    u, (q, p)
end


struct Bind{T,S} <: DynamicIterator
    Y::S
    P::T
end

function dyniterate(M::Bind, start::Start)
    v, q = @returnnothing iterate(M.Y)
    u, p = @returnnothing dyniterate(M.P, start, v)
    u, (q, p)
end

function iterate(M::Bind, (q, p)::Tuple)
    v, q = @returnnothing iterate(M.Y, q)
    u, p = @returnnothing dyniterate(M.P, p, v)
    u, (q, p)
end


struct ControlledFilter{R,T,S} <: DynamicIterator
    C::R
    Y::S
    P::T
end
filter(Y, P; control = nothing) = filter(Y, P, control)
filter(Y, P, C) = ControlledFilter(C, Y, P)



function dyniterate(M::ControlledFilter, ::Nothing, (value,)::Value)
    w, r = @returnnothing dyniterate(M.C, nothing, (control = value,))
    v, q = @returnnothing iterate(M.Y)
    u, p = @returnnothing dyniterate(M.P, nothing, (value = value, control = w, observation = v))
    u, (u, r, q, p)
end

function iterate(M::ControlledFilter, (u, r, q, p)::Tuple)
    w, r = @returnnothing iterate(M.C, r, (control = u,))
    v, q = @returnnothing iterate(M.Y, q)
    u, p = @returnnothing dyniterate(M.P, p, (control = w, observation = v))
    u, (u, r, q, p)
end
