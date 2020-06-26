var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Kalman.jl-1",
    "page": "Home",
    "title": "Kalman.jl",
    "category": "section",
    "text": "A Julia package for fast, flexible filtering and smoothing.At the moment implements the vanilla Kalman filter, but the layout allows for easy extension. I took some effort to allow for interoperability with StaticArrays and ForwardDiff.See the Library tab for the assorted documentation."
},

{
    "location": "library/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library/#Library-1",
    "page": "Library",
    "title": "Library",
    "category": "section",
    "text": ""
},

{
    "location": "library/#Kalman.LinearObservation",
    "page": "Library",
    "title": "Kalman.LinearObservation",
    "category": "type",
    "text": "LinearObservation(P, H, R)\n\nObserve the LinearEvolution P using y = Hx + v where v sim N(0 R).\n\nExamples\n\n    O = LinearObservation(LinearEvolution(Φ, Gaussian(b, Q)), H, R)\n\n\n\n\n\n"
},

{
    "location": "library/#Kalman.LinearEvolution",
    "page": "Library",
    "title": "Kalman.LinearEvolution",
    "category": "type",
    "text": "LinearEvolution(Φ, b, Q) <: Evolution\n\nEvolution of the law of x -> Φ x + w where w sim N(0 Q).\n\nExamples\n\n    evolve(LinearEvolution(Φ, b, Q), 0 => Gaussian(x, P))\n\n\n\n\n\n"
},

{
    "location": "library/#Kalman.LinearStateSpaceModel",
    "page": "Library",
    "title": "Kalman.LinearStateSpaceModel",
    "category": "type",
    "text": "LinearStateSpaceModel <: StateSpaceModel\n\nLinearStateSpaceModel(sys, obs)\n\nCombines a linear system sys and an observations model obs and to a linear statespace model in a modular way.\n\nEvolves StateObs objects.\n\n\n\n\n\n"
},

{
    "location": "library/#Linear-systems-1",
    "page": "Library",
    "title": "Linear systems",
    "category": "section",
    "text": "Kalman.LinearObservation\nKalman.LinearEvolution\nKalman.LinearStateSpaceModel"
},

{
    "location": "library/#Kalman.kalmanfilter",
    "page": "Library",
    "title": "Kalman.kalmanfilter",
    "category": "function",
    "text": "kalmanfilter(M, t => x0) -> kf\n\nkf(iter) is an iterator over Gaussians or Distributions representing the filtered distribution of x where y iterates over (enumerated) signal values.\n\nExample\n\nkf = kalmanfilter(M, 0 => prior) #\nest1 = collect(kf(Y1))\nest2 = collect(kf(Y2))\n\n\n\n\n\n"
},

{
    "location": "library/#Kalman.Filtered",
    "page": "Library",
    "title": "Kalman.Filtered",
    "category": "type",
    "text": "filter(Y, P)\n\n\"Filter\" data Y with iterator P calling (handling of nothings omitted)\n\n(t, y), state = iterate(Y, state)\n(s => x) = dyniterate(P, s => x, (observation = t => y,))\n\n\n\n\n\n"
},

{
    "location": "library/#Filtering-1",
    "page": "Library",
    "title": "Filtering",
    "category": "section",
    "text": "Kalman.kalmanfilter\nKalman.Filtered"
},

{
    "location": "library/#Smoothing-1",
    "page": "Library",
    "title": "Smoothing",
    "category": "section",
    "text": ""
},

{
    "location": "library/#Iterators-1",
    "page": "Library",
    "title": "Iterators",
    "category": "section",
    "text": "DynamicIterators.dyniterate"
},

]}
