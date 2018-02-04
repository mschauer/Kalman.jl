
"""
    xxhat, Mhat = kalmanEM(yy, M1::LinearHomogSystem, maxit = 500, tol = 0.001) 

EM procedure for the state space model (Shumway, Stoffer)

    yy -- data
    M1 -- initial model
    maxit -- maximum number of iterations
    tol -- rel. tolerance for test for convergence of marginal likelihood


    xxhat, Mhat -- smoothed process with estimated model
"""
function kalmanEM(yy, M::LinearHomogSystem, maxit = 500, tol = 0.001) 

    H, Phi, b, Q, R = M.H, M.Phi, M.b, M.Q, M.R
    d2, d = size(M.H)
    assert(ndims(yy) == 2)
    assert(norm(b) == 0)
    n = size(yy,2)

    x0 = M.x0
    P0 = M.P0
    xxs = xxf = zeros(d, n) # #x(f/s)[k] predited and corrected (and smoothed)
       
    PPf = zeros(d, d, n)  #P[k,k]
    PPpred = zeros(d, d, n)     #P[k, k-1]
    PPs = zeros(d, d, n)  #P[k,k]
    PPcs = zeros(d, d, n)  #P[k-1,k]

        
    K = zeros(d, d)
    
    cvg = 1.+tol
    nll = NaN
    for k in 1:maxit
        xf = x0
        Pf = P0

        # E step, forward pass
        nllold = nll
        nll = 0.
        for i in 1:n
            xf, Pf, Ppred, l, K  = kalman_kernel(xf, yy[:, i], Pf, H, Phi, b, Q, R)
            xxf[:, i], PPf[:, :, i], PPpred[:, :, i] = xf, Pf, Ppred  
            l += d*log(2pi)/2         # constant messes up convergence test      
            nll -= l
        end
    
        if k > 1
            cvg = (nllold-nll)/abs(nllold)
            cvg < 0 && warn("Likelihood decreasing")
            if abs(cvg) < tol
#                println([nll, cvg])
                break
            end
        end
        
        # E step, backwards pass including PPcs

        #start with xf, Pf from forward pass
        xs = xf
        Ps = Pf

        PPcs[:, :, n] = (I - K*H)*Phi*PPf[:, :, n-1]
        J = PPf[:, :, n]*Phi'/(Phi*Ps*Phi' + Q)
        PPs[:, :, n] = Ps
        for i in n-1:-1:1
            Jprev = J   
            xs, Ps, J = smoother_kernel(xs, Ps, xxf[:, i],  PPf[:, :, i], PPpred[:, :, i+1], Phi, b)
            xxs[:, i], PPs[:, :, i] = xs, Ps
            
            j = i + 2 #i = j - 2 as J = J(i)
            if 2 < j <= n
                
                 Jjm2 = J
                 Jjm1 = Jprev
                 PPcs[:, :, j-1] = PPf[:, :, j-1]*Jjm2' +  Jjm1*(PPcs[:, :, j]  - Phi*PPf[:, :, j-1])*Jjm2'
                 #Pcs[,,j-1]=Pf[,,j-1]%*%t(J[,,j-2])+ J[,,j-1]%*%(Pcs[,,j]-Phi%*%Pf[,,j-1])%*%t(J[,,j-2])}
                 
            end
        end
        x0, P0, J0 = smoother_kernel(xs, Ps, x0, P0, PPpred[:, :, 1], Phi, b)
        PPcs[:, :, 1] = PPf[:, :, 1]*J0' +  J*(PPcs[:, :, 2]  - Phi*PPf[:, :, 1])*J0'
    

        # M step

        A11 = xxs[:, 1]*xxs[:, 1]' + PPs[:, :, 1]
        A10 = xxs[:, 1]*x0' + PPcs[:, :, 1] 
        A00 = x0*x0' + P0
        
        u = yy[:,1]-H*xxs[:,1]
        R = u*u' + H*PPs[:, :, 1]*H'
        
        for i in 1:n-1
            A11 = A11 + xxs[:, i+1]*xxs[:, i+1]' + PPs[:, :, i+1]
            A10 = A10 + xxs[:, i+1]*xxs[:, i]' + PPcs[:, :, i+1] 
            A00 = A00 + xxs[:, i]*xxs[:, i]' + PPs[:, :, i]
            
            u = yy[:,i+1]-H*xxs[:,i+1]
            R = R + u * u' + H*PPs[:, :, i+1]*H'
        end
      
        Phi = A10/A00
        Q = (A11 - Phi*A10')/n
        Q = (Q+Q')/2
        R = R/n

     
#        println("PARAM $k: nll", [nll, cvg], "\n x0 $x0\n P0 $P0\n R $R\n Phi $Phi\n Q $Q")

    end

    
    xxs, LinearHomogSystem(x0, P0, Phi, b, Q, H, R)
end
