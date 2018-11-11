
struct GenericLinearObservation <: AbstractObservation
end


"""
    LinearObservation(H, R)

Observe `y = Hx + v` where ``v ~ N(0, R)``
"""
struct LinearObservation{TH, TR} <: AbstractObservation
    H::TH # dxd
    R::TR # d
end
