
x0 = [1., 0.]
P0 = eye(2)

Phi = [0.8 0.5; 0.0 0.8]
b = zeros(2)
Q = [0.1 0.0; 0.0 1.0]

yshadow = [1.0]
H = [1.0 0.0]
R = eye(1)
M = LinearHomogSystem(x0, P0, Phi, b, Q, yshadow, H, R)
