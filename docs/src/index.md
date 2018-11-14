# Kalman.jl

A Julia package for fast, flexible filtering and smoothing.

At the moment implements the vanilla Kalman filter, but the layout allows for easy extension. I took some effort to allow for interoperability with `StaticArrays` and `ForwardDiff`.

```@docs
Kalman.LinearEvolution
Kalman.LinearObservationModel
Kalman.kalmanfilter!
```@docs

See the Library tab for the assorted documentation.
