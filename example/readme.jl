# Define linear evolution
Φ = [0.8 0.5; -0.1 0.8]
b = zeros(2)
Q = [0.2 0.0; 0.0 1.0]

E = LinearEvolution(Φ, Gaussian(b, Q))

# Define observation scheme
H = [1.0 0.0]
R = Matrix(1.0I, 1, 1)

O = LinearObservation(E, H, R)

# Prior
x0 = [1., 0.]
P0 = Matrix(1.0I, 2, 2)
prior = 0 => Gaussian(x0, P0)

# Observations (mock)
Y = [1 => [1.14326], 2 => [-0.271804], 3 => [-0.00512675]]

# Filter
Xf, ll = kalmanfilter(O, prior, Y)
@show Xf



# Online
@testset "Online" begin
    ϕ = iterate(Y)
    ϕ === nothing && error("no observations")
    y, ystate = ϕ

    ϕ = dyniterate(O, Start(Kalman.Filter(prior, 0.0)), y)
    ϕ === nothing && error("no observations")
    (t, u), state = ϕ

    X = trajectory((t => u[1],))
    while true
        ϕ = iterate(Y, ystate)
        ϕ === nothing && break
        y, ystate = ϕ

        ϕ = dyniterate(O, state, y)
        ϕ === nothing && break
        (t, u), state = ϕ
        push!(X, t => u[1]) # filtered state as Gaussian
    end
    ll = u[3] # likelihood
    @show  X, ll
end
