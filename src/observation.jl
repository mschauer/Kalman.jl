

abstract type Observation
end

struct GenericLinearObservation <: Observation
end

"""
    LinearObservation(H, R)

Observe `y = Hx + v` where ``v ~ N(0, R)``
"""
struct LinearObservation{TH, TR} <: Observation
    H::TH # dxd
    R::TR # d
end
