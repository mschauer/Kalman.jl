var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Kalman.jl-1",
    "page": "Home",
    "title": "Kalman.jl",
    "category": "section",
    "text": "A Julia package for fast, flexible filtering and smoothing."
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Library-1",
    "page": "Library",
    "title": "Library",
    "category": "section",
    "text": ""
},

{
    "location": "library.html#Kalman.LinearHomogSystem",
    "page": "Library",
    "title": "Kalman.LinearHomogSystem",
    "category": "Type",
    "text": "LinearHomogSystem\n\nLinear dynamical model with linear observation scheme\n\n\n\n"
},

{
    "location": "library.html#StatsBase.sample",
    "page": "Library",
    "title": "StatsBase.sample",
    "category": "Function",
    "text": "Y, X = sample(n, m, M::LinearHomogSystem) \n\nSample observations\n\nn – number of observations per simulation m – number of independent simulations M – linear dynamic system and observation model\n\nY – observations (d2xnxm array) X – simulated processes (dxnxm array)\n\n\n\n"
},

{
    "location": "library.html#Linear-systems-1",
    "page": "Library",
    "title": "Linear systems",
    "category": "section",
    "text": "Kalman.LinearHomogSystem\nStatsBase.sample"
},

{
    "location": "library.html#Kalman.kalmanfilter",
    "page": "Library",
    "title": "Kalman.kalmanfilter",
    "category": "Function",
    "text": "ll, xxf = kalmanfilter(yy::Matrix, M::LinearHomogSystem)\n\nKalman filter\n\nyy -- d₂xn array\nM -- linear model\nll -- marginal likelihood\nxxf -- filtered process\n\n\n\nll, X = kalmanfilter{T}(Y::Array{T,3}, M::LinearHomogSystem)\n\nStack Kalman filter\n\nY -- array of m independent processes (dxnxm array)\n\nX -- filtered processes (dxnxm array)\nll -- product marginal log likelihood\n\n\n\n"
},

{
    "location": "library.html#Kalman.observe!",
    "page": "Library",
    "title": "Kalman.observe!",
    "category": "Function",
    "text": "observe!(s, x, P, t, y, M::LinearHomogSystem) -> t, y, H, R\n\n\n\n"
},

{
    "location": "library.html#Kalman.predict!",
    "page": "Library",
    "title": "Kalman.predict!",
    "category": "Function",
    "text": "predict!(s, x, P, t, M::LinearHomogSystem) -> x, Ppred, Phi\n\n\n\n"
},

{
    "location": "library.html#Kalman.correct!",
    "page": "Library",
    "title": "Kalman.correct!",
    "category": "Function",
    "text": "correct!(x, Ppred, y, H, R, SSM) -> x, P, yres, S, K\n\n\n\n"
},

{
    "location": "library.html#Kalman.evolve",
    "page": "Library",
    "title": "Kalman.evolve",
    "category": "Function",
    "text": "evolve(s, x, P, t, M::LinearHomogSystem) -> x\n\n\n\n"
},

{
    "location": "library.html#Kalman.kalman_kernel",
    "page": "Library",
    "title": "Kalman.kalman_kernel",
    "category": "Function",
    "text": "kalman_kernel(s, x, P, t, Y, SSM) -> t, x, P, Ppred, ll, K\n\nSingle Kalman filter step.\n\n\n\n"
},

{
    "location": "library.html#Filtering-1",
    "page": "Library",
    "title": "Filtering",
    "category": "section",
    "text": "Kalman.kalmanfilter\nKalman.observe!\nKalman.predict!\nKalman.correct!\nKalman.evolve\nKalman.kalman_kernel"
},

{
    "location": "library.html#Kalman.kalmanrts",
    "page": "Library",
    "title": "Kalman.kalmanrts",
    "category": "Function",
    "text": "xxs, PP, PPpred, ll = kalmanrts(yy, xxs, M::LinearHomogSystem)\n\nRauch-Tung-Striebel smoother\n\nyy -- d₂xn array or \nxxs -- empty dxn array or view\nM --  linear dynamic system and observation model\n\nxxs -- smoothed process\nPP -- smoother variance\nPPpred -- smoother prediction\nll -- filter likelihood\n\n\n\nX = kalmanrts{T}(Y::Array{T,3}, M::LinearHomogSystem)\n\nStack Rauch-Tung-Striebel smoother, computes the marginal smoothed states, that is it computes the law p(x_i  y_1n).\n\nY -- m independent processes (d₂xnxm array)\n\nX -- smoothed processes (dxnxm array)\n\n\n\n"
},

{
    "location": "library.html#Smoothing-1",
    "page": "Library",
    "title": "Smoothing",
    "category": "section",
    "text": "Kalman.kalmanrts"
},

{
    "location": "library.html#Kalman.KalmanFilter",
    "page": "Library",
    "title": "Kalman.KalmanFilter",
    "category": "Type",
    "text": "KalmanFilter(y, M)\n\nKalman filter as iterator, iterating over Gaussians representing the filtered distribution of x. Arguments y iterates over signal values.\n\n\n\n"
},

{
    "location": "library.html#Kalman.TimedKalmanFilter",
    "page": "Library",
    "title": "Kalman.TimedKalmanFilter",
    "category": "Type",
    "text": "TimedKalmanFilter(ty, M)\n\nKalman filter as iterator, iterating over Gaussians representing the filtered distribution of x. ty iterates over pairs (t, y)  of signal time and signal value.\n\n\n\n"
},

{
    "location": "library.html#Iterators-1",
    "page": "Library",
    "title": "Iterators",
    "category": "section",
    "text": "Kalman.KalmanFilter\nKalman.TimedKalmanFilter"
},

{
    "location": "library.html#Kalman.track",
    "page": "Library",
    "title": "Kalman.track",
    "category": "Function",
    "text": "track(A::Matrix{Float64}, p, h = 10) -> p2, err\n\nTrack blurry lightsource by applying a window with half-width h at  an approximate location p = (i,j) and find the average weighted location of points with high light intensity. \n\nGives an standard error estimate.\n\n\n\n"
},

{
    "location": "library.html#Tracking-1",
    "page": "Library",
    "title": "Tracking",
    "category": "section",
    "text": "Kalman.track"
},

]}
