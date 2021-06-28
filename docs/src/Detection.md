# Detection

## ECG

```@docs
detectRPeaks
```

```@setup rpeaks
# Random ECG of length 40s with a known number of 52 Beats
using Plots, Cardio, DataFrames, CSV
ecg = CSV.read("../data/ecg.csv", DataFrame)[!, :ECG]
```

```@example rpeaks
plot(ecg, lab = "")
peaks = detectRPeaks(ecg, 250) # signal was sampled at 250 Hz
scatter!(peaks, ecg[peaks])
```

```@docs
getECGBaseline
```

```@setup baseline
# Random ECG of length 40s with a known number of 52 Beats
using Plots, Cardio, DataFrames, CSV
ecg = CSV.read("../data/ecg.csv", DataFrame)[!, :ECG]
```

```@example baseline
plot(ecg, lab = "")
baseline = getECGBaseline(ecg, 250) # signal was sampled at 250 Hz
plot!(baseline, lab = "baseline")
```
## Blood Pressure