
x0 = 1.0
P0 = 1.0

Φ = 0.5
b = 0.0
Q = 2.0

yshadow = 0.0
H = 1.0
R = 1.0

E = LinearEvolution(Φ, Gaussian(b, Q))
Obs = LinearObservationModel(H, R)
M = LinearStateSpaceModel(E, Obs)

O = LinearObservation(LinearEvolution(Φ, Gaussian(b, Q)), H, R)
