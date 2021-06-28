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

