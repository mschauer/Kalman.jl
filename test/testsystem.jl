
x0 = 1.0
P0 = 1.0

Φ = -1.0
b = 0.0
Q = 1.0

yshadow = 0.0
H = 1.0
R = 1.0

E = LinearEvolution(Φ, Gaussian(b, Q))
Obs = LinearObservation(H, R)
M = LinearStateSpaceModel(E, Obs)
