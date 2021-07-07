# BRS
Several methods to estimate baroreflex sensitivity (BRS) have been defined, many of which are contained in this package. All methods require the following inputs:

* RR: The RR Interval series in ms
* SBP: The systolic blood pressure in mmHg

Several methods can be tweaked by using keywords, though the most common and recommended parameters are set as default.

```@docs
BRS.getBRS
```

## Sequence Method

```@docs
BRS.sme
```

```@docs
BRS.SME
```

```@setup sme
using Plots, Cardio, DataFrames, CSV
input = CSV.read("../data/BRS.csv", DataFrame)
```

```@example sme
result = BRS.sme(input.RR, input.SBP)
plot(result, dpi = 120)
```

## xBRS Method
```@docs
BRS.xbrs
```

```@docs
BRS.xBRS
```

```@setup xBRS
using Plots, Cardio, DataFrames, CSV
input = CSV.read("../data/BRS.csv", DataFrame)
```

```@example xBRS
result = BRS.xbrs(input.RR, input.SBP)
plot(result, dpi = 120)
```

## RMSSD Ratio
```@docs
BRS.rmssdr
```

## Transfer Function Method

```@docs
BRS.tfbrs
```

```@docs
BRS.tfBRS
```

```@setup tfBRS
using Plots, Cardio, DataFrames, CSV
input = CSV.read("../data/BRS.csv", DataFrame)
```

```@example tfBRS
result = BRS.tfbrs(input.RR, input.SBP, n = 100)
plot(result, dpi = 120)
```

## Spectral AR-based Method

```@docs
BRS.arbrs
```

```@docs
BRS.arBRS
```

```@setup arBRS
using Plots, Cardio, DataFrames, CSV, Statistics
input = CSV.read("../data/BRS.csv", DataFrame)
```

```@example arBRS
result = BRS.arbrs(input.RR, input.SBP; sf = 1000 / mean(input.RR))
plot(result, dpi = 120)
```

## Phase-Rectified Signal Averaging Method
Source: Bauer, A., Morley-Davies, A., Barthel, P., Müller, A., Ulm, K., Malik, M., & Schmidt, G. (2010). Bivariate phase-rectified signal averaging for assessment of spontaneous baroreflex sensitivity: pilot study of the technology. Journal of electrocardiology, 43(6), 649-653. [DOI](https://doi.org/10.1016/j.jelectrocard.2010.05.012)
```@docs
BRS.prsabrs
```

```@docs
BRS.prsaBRS
```

```@setup prsabrs
using Plots, Cardio, DataFrames, CSV
input = CSV.read("../data/BRS.csv", DataFrame)
```

```@example prsabrs
result = BRS.prsabrs(input.RR, input.SBP)
plot(result, dpi = 120)
```