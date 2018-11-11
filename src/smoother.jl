
function smoother_kernel(Gs::T, Gf, Ppred, Phi, b) where {T}
    xs, Ps = meancov(Gs)
    xf, Pf = meancov(Gf)

    J = Pf*Phi'/Ppred # C/(C+w)
    xs = xf +  J*(xs - (Phi*xf  + b))
    Ps = Pf + J*(Ps - Ppred)*J'

    T(xs, Ps), J
end
