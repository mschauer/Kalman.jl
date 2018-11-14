x0 = [1., 0.]
P0 = Matrix(1.0I, 2, 2)

Φ = [0.8 0.5; -0.1 0.8]
b = zeros(2)
Q = [0.2 0.0; 0.0 1.0]

yshadow = [0.0]
H = [1.0 0.0]
R = Matrix(1.0I, 1, 1)

E = LinearEvolution(Φ, Gaussian(b, Q))
Obs = LinearObservationModel(H, R)
M = LinearStateSpaceModel(E, Obs)
