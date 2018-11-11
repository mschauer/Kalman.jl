
function backward_kernel(Gs::T, Gf, Ppred, Phi, b) where T
    xs, Ps = meancov(Gs)
    xf, Pf = meancov(Gf)

    J = Pf*Phi'/Ppred
    xs = xf +  J*(x - (Phi*xf  + b))
    Ps = Pf - J*Ppred*J' # as Ppred = M.Phi*P*M.Phi' + M.Q
    T(xs, Ps), J
end
