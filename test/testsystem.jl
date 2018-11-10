
x0 = 1.0
P0 = 1.0

Phi = -1.0
b = 0.0
Q = 1.0

yshadow = 0.0
H = 1.0
R = 1.0
M = LinearHomogSystem(Gaussian(x0, P0), Phi, Gaussian(b, Q), H, Gaussian(yshadow, R))
