var documenterSearchIndex = {"docs":
[{"location":"HRV/#HRV","page":"HRV","title":"HRV","text":"","category":"section"},{"location":"Utility/#Utility-Functionality","page":"Utility","title":"Utility Functionality","text":"","category":"section"},{"location":"Utility/","page":"Utility","title":"Utility","text":"adaptiveHRVFilter","category":"page"},{"location":"Utility/#Cardio.adaptiveHRVFilter","page":"Utility","title":"Cardio.adaptiveHRVFilter","text":"adaptiveHRVFilter(signal::Vector{<:Real}; removeoutliers::Bool = true, replacenonnormal::Bool = true, adaptivecontrollingcoef::Real = 0.05, proportionalitylimit::Real = 10/100, outlierminzfactor::Real = 3, maxexcesshrv::Real = 20, physiologicalvalues::Tuple{Real, Real} = (200, 2000))\n\nAn addaptive filter for HRV data      Based on: Wessel, N., Voss, A., Malberg, H., Ziehmann, Ch.,          Voss, H. U., Schirdewan, A., Meyerfeldt, U.,Kurths, J.:          Nonlinear analysis of complex phenomena in cardiological data,          Herzschr. Elektrophys., 11(3), 2000, 159-173, doi:10.1007/s003990070035.\n\nFilters out unphysiological beats in the RR series with the ability to replace them.\n\nArgs:\n\n'signal::Vector': HRV in ms\n\nKeyword Args\n\n'remove_outliers::Bool = true': option if non physiological outliers shall be removed\n'replace_nonnormal::Bool = true': if true, non normal HRV values are replaced\n'adaptivecontrollingcoef::Real = 0.05': ???\n'proportionality_limit::Real = 10/100': ???\n'outlierminzfactor::Real = 3': ???\n'maxexcesshrv::Real = 20': ???\n'physiological_values::Tuple{Real, Real} = (200, 2000)': Definition of the physiological values, scheme: (min, max)\n\nReturn\n\nReturn value depends on keyword args\n\njulia> adaptiveHRVFilter(signal)\n\n\n\n\n\n","category":"function"},{"location":"Utility/","page":"Utility","title":"Utility","text":"medfilt1","category":"page"},{"location":"Utility/#Cardio.medfilt1","page":"Utility","title":"Cardio.medfilt1","text":"medfilt1(x::Array{<:Real}; n::Int = 3, padding::String = \"zeropad\", dim::Int = -1)\n\nApply a median filter to a signal vector or array x, similar to Matlabs medfilt1. Using Heap based calculation of the median to increase performance for larger windows n.\n\nArgs:\n\n'x::Array{<:Real}': Array containing real values\n\nKeywords:\n\n'n::Int': Window length. The Median at point i is defined as median(x[i-n+1:i])\n'padding::String': Specifies how to deal with Endpoints. The modes 'zeropad' and 'truncate' are available, with the first as default.\n'dim::Int': Specifies the dimension to be filtered along. As default the first non singleton dimension is chosen.\n\nReturn:\n\n'Array{Float64,N}': Always type Float64 with the same length as the input x\n\nExamples\n\njulia> medfilt1(collect(1:10))\n10-element Array{Float64,1}:\n 1.0\n 2.0\n 3.0\n 4.0\n 5.0\n 6.0\n 7.0\n 8.0\n 9.0\n 9.0\n\n\n\n\n\n","category":"function"},{"location":"BRS/#BRS","page":"BRS","title":"BRS","text":"","category":"section"},{"location":"BRS/","page":"BRS","title":"BRS","text":"Several methods to estimate baroreflex sensitivity (BRS) have been defined, many of which are contained in this package. All methods require the following inputs:","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"RR: The RR Interval series in ms\nSBP: The systolic blood pressure in mmHg","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"Several methods can be tweaked by using keywords, though the most common and recommended parameters are set as default.","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.getBRS","category":"page"},{"location":"BRS/#Cardio.BRS.getBRS","page":"BRS","title":"Cardio.BRS.getBRS","text":"    getBRS(RR::Vector{<:Real}, SBP::Vector{<:Real})\n\nReturns all BRS measures with default settings for individual methods.\n\nMethods\n\nSME: Sequence Method\nRMSSDR: RMSSD ratio\nxBRS: Cross-correlation baroreflex sensitivity\n\n\n\n\n\n","category":"function"},{"location":"BRS/#Sequence-Method","page":"BRS","title":"Sequence Method","text":"","category":"section"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.sme","category":"page"},{"location":"BRS/#Cardio.BRS.sme","page":"BRS","title":"Cardio.BRS.sme","text":"sme(RR::Vector{<:Real}, SBP::Vector{<:Real}; thresholdRR::Float64 = 4.0, thresholdSBP::Float64 = 1.0, seqLen::Int = 3, delay::Int = 1, minCor::Float64 = 0.8)\n\nCalculate the BaroReflex Sensitivity (BRS) for a serious of RR intervals and respectiv systolic bloodpressure using the Sequence Method (SME).\n\nArgs:\n\n'RR::Vector': Data Vector containing the RR intervals\n'SBP::Vector': Data Vector containing the systolic bloodpressure\n\nKeywords:\n\n'thresholdRR::Float64': Threshold for change in RR interval to qualify for a valid sequence (literature: 4ms)\n'thresholdSBP::Float64': Threshold for change in SBP to qualify for a valid sequence (literature: 1mmHg)\n'seqLen::Int': Minimum length of a valid sequence (literature: 3)\n'delay::Int': Delay between RR and SBP (literature: 1)\n'minCor::Float64': The minimal correlation between RR and SBP in a sequence to qualify as a valid sequence\n'BRS::Float64': The BRS as estimated by the Sequence Method\n\nReturn:\n\nReturns a SME struct for which a plotting recipe is provided, so that found sequences can be inspected by calling plot() when using Plots.jl\n\n\n\n\n\n","category":"function"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.SME","category":"page"},{"location":"BRS/#Cardio.BRS.SME","page":"BRS","title":"Cardio.BRS.SME","text":"Struct that stores all information related to the sequence method. The main result is stored in 'sBRS'. It can be plotted for visual inspection.\n\n\n\n\n\n","category":"type"},{"location":"BRS/","page":"BRS","title":"BRS","text":"using Plots, Cardio, DataFrames, CSV\ninput = CSV.read(\"../data/BRS.csv\", DataFrame)","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"result = BRS.sme(input.RR, input.SBP)\nplot(result, dpi = 120)","category":"page"},{"location":"BRS/#xBRS-Method","page":"BRS","title":"xBRS Method","text":"","category":"section"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.xbrs","category":"page"},{"location":"BRS/#Cardio.BRS.xbrs","page":"BRS","title":"Cardio.BRS.xbrs","text":"xbrs(RR::Vector{<:Real}, SBP::Vector{<:Real};...)\n\nCalculate the xBRS index for the assesment of baroreflex sensitivity. Based on Westerhof, B. E. et al. (2004). Time-domain cross-correlation baroreflex sensitivity: performance on the EUROBAVAR data set. Journal of hypertension, 22(7), 1371-1380\n\nKeyword Arguments\n\nminCor: Minimal significant correlation, dafaults to 0.632 (p = 0.05, two-sided for 10 s window)\ntExcerpt: Length of the sliding window in seconds, defaults to 10 s\ndelays: Which delays, shift of the RR window shall be considered, defaults to 0:5\n\nReturn\n\nReturns an xBRS structure, that can be be plotted for visual inspection.\n\n\n\n\n\n","category":"function"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.xBRS","category":"page"},{"location":"BRS/#Cardio.BRS.xBRS","page":"BRS","title":"Cardio.BRS.xBRS","text":"Struct that stores all information regarting the xBRS etimation. The final result is stored in 'xBRSg'. It can be plotted for visual inspection.\n\n\n\n\n\n","category":"type"},{"location":"BRS/","page":"BRS","title":"BRS","text":"using Plots, Cardio, DataFrames, CSV\ninput = CSV.read(\"../data/BRS.csv\", DataFrame)","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"result = BRS.xbrs(input.RR, input.SBP)\nplot(result, dpi = 120)","category":"page"},{"location":"BRS/#RMSSD-Ratio","page":"BRS","title":"RMSSD Ratio","text":"","category":"section"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.rmssdr","category":"page"},{"location":"BRS/#Cardio.BRS.rmssdr","page":"BRS","title":"Cardio.BRS.rmssdr","text":"rmssdr(RR::Vector{<:Real}, SBP::Vector{<:Real})\n\nReturns the RMSSD ratio := RMSSD(RR) / RMSSD(SBP)\n\n\n\n\n\n","category":"function"},{"location":"BRS/#Transfer-Function-Method","page":"BRS","title":"Transfer Function Method","text":"","category":"section"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.tfbrs","category":"page"},{"location":"BRS/#Cardio.BRS.tfbrs","page":"BRS","title":"Cardio.BRS.tfbrs","text":"tfbrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; ...)\n\nTransfer function based BRS measure as defined by Robbe et al. \n\nKeyword Arguments\n\nn: length of hamming window for spectral estimation, defaults to length(RR) ÷ 10\nminCoh: minimal valid coherence, defaults to 0.5\nLF: The frequency range defined as low frequeny, defaults to (0.04, 0.15)\n\n\n\n\n\n","category":"function"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.tfBRS","category":"page"},{"location":"BRS/#Cardio.BRS.tfBRS","page":"BRS","title":"Cardio.BRS.tfBRS","text":"Struct that stores all information regarting the tfBRS etimation. The final result is stored in 'tfBRSv'. It can be plotted for visual inspection.\n\n\n\n\n\n","category":"type"},{"location":"BRS/","page":"BRS","title":"BRS","text":"using Plots, Cardio, DataFrames, CSV\ninput = CSV.read(\"../data/BRS.csv\", DataFrame)","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"result = BRS.tfbrs(input.RR, input.SBP, n = 100)\nplot(result, dpi = 120)","category":"page"},{"location":"BRS/#Phase-Rectified-Signal-Averaging-Method","page":"BRS","title":"Phase-Rectified Signal Averaging Method","text":"","category":"section"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.prsabrs","category":"page"},{"location":"BRS/#Cardio.BRS.prsabrs","page":"BRS","title":"Cardio.BRS.prsabrs","text":"prsabrs(RR::Vector{<:Real}, SBP::Vector{<:Real}; L::Int = 15)\n\nCalculate a measure for BRS based on phase-rectified signal averaging as defined by Bauer et al. 2010.\n\nKeywords\n\nL: defines the segment length as 2L+1\n\n\n\n\n\n","category":"function"},{"location":"BRS/","page":"BRS","title":"BRS","text":"BRS.prsaBRS","category":"page"},{"location":"BRS/#Cardio.BRS.prsaBRS","page":"BRS","title":"Cardio.BRS.prsaBRS","text":"Struct that stores all information regarting the prsaBRS etimation. The final result is stored in 'prsaBRSv'. It can be plotted for visual inspection.\n\n\n\n\n\n","category":"type"},{"location":"BRS/","page":"BRS","title":"BRS","text":"using Plots, Cardio, DataFrames, CSV\ninput = CSV.read(\"../data/BRS.csv\", DataFrame)","category":"page"},{"location":"BRS/","page":"BRS","title":"BRS","text":"result = BRS.prsabrs(input.RR, input.SBP, n = 100)\nplot(result, dpi = 120)","category":"page"},{"location":"Detection/#Detection","page":"Detection","title":"Detection","text":"","category":"section"},{"location":"Detection/#ECG","page":"Detection","title":"ECG","text":"","category":"section"},{"location":"Detection/","page":"Detection","title":"Detection","text":"detectRPeaks","category":"page"},{"location":"Detection/#Cardio.detectRPeaks","page":"Detection","title":"Cardio.detectRPeaks","text":"detectRPeaks(ecg::Vector{<:Real}, samplerate::Real; minPeakDist::Real = 0.360)\n\nFind R peaks in ECG signals as specified by Benitez et al. See http://dx.doi.org/10.1016/S0010-4825(01)00009-9 for more information.\n\nArgs:\n\necg: ECG data \nsamplerate: Sampling rate [Hz]\nminPeakDist: minimum distance between consecutive peaks [s]\n\nReturn:\n\n'res::Vector{Int64}': Vector containing the position of the R peaks in ecg, divide by samplerate to get values in a time base\n\n\n\n\n\n","category":"function"},{"location":"Detection/","page":"Detection","title":"Detection","text":"# Random ECG of length 40s with a known number of 52 Beats\nusing Plots, Cardio, DataFrames, CSV\necg = CSV.read(\"../data/ecg.csv\", DataFrame)[!, :ECG]","category":"page"},{"location":"Detection/","page":"Detection","title":"Detection","text":"plot(ecg, lab = \"\")\npeaks = detectRPeaks(ecg, 250) # signal was sampled at 250 Hz\nscatter!(peaks, ecg[peaks], lab = \"R peaks\")","category":"page"},{"location":"Detection/","page":"Detection","title":"Detection","text":"getECGBaseline","category":"page"},{"location":"Detection/#Cardio.getECGBaseline","page":"Detection","title":"Cardio.getECGBaseline","text":"getECGBaseline(ecg::Vector{<:Real}, samplerate::Real)\n\nGet the baseline of an ECG signal for baseline correction. Source: Advances in Cardiac Signal Processing - Acharya, U.R. and Suri, J. and Spaan, J.A.E. and Krishnan, S.M. and Technologies, B.      - ISBN: 9783540366751 page: 58f. adaption by Jan F. Kraemer\n\nArgs:\n\n'ecg::Vector{<:Number}': ECG signal\n'samplerate::Number': Sampling rate\n\n\n\n\n\n","category":"function"},{"location":"Detection/","page":"Detection","title":"Detection","text":"# Random ECG of length 40s with a known number of 52 Beats\nusing Plots, Cardio, DataFrames, CSV\necg = CSV.read(\"../data/ecg.csv\", DataFrame)[!, :ECG]","category":"page"},{"location":"Detection/","page":"Detection","title":"Detection","text":"plot(ecg, lab = \"\")\nbaseline = getECGBaseline(ecg, 250) # signal was sampled at 250 Hz\nplot!(baseline, lab = \"baseline\", linewidth = 2)","category":"page"},{"location":"Detection/#Blood-Pressure","page":"Detection","title":"Blood Pressure","text":"","category":"section"},{"location":"Detection/","page":"Detection","title":"Detection","text":"detectPWPeaks","category":"page"},{"location":"Detection/#Cardio.detectPWPeaks","page":"Detection","title":"Cardio.detectPWPeaks","text":"detectPWPeaks(signal::Vector{<:Real}, fs::Real;...)\n\nFind Peaks in a Pulswave signal with an algorithm proposed by Nenova, B., & Iliev, I. (2010). An automated algorithm for fast pulse wave detection. International Journal Bioautomation, 14(3), 203.\n\nArgs:\n\nsignal: Pulswave data\nfs: Sampling Frequency [Hz]\nwindowLenght: Optional lenght of the window in ms. (standard: 300ms)\ntuning: fine tuning to filter Peaks below a certain height\n\nReturn:\n\nReturns a Vector containing the Peak indices. Use signal[indices] to access the Peak values.\n\n\n\n\n\n","category":"function"},{"location":"Detection/","page":"Detection","title":"Detection","text":"using Plots, Cardio, DataFrames, CSV\nbp = CSV.read(\"../data/BP.csv\", DataFrame)[!, :BP]","category":"page"},{"location":"Detection/","page":"Detection","title":"Detection","text":"plot(bp, lab = \"\")\npeaks = detectPWPeaks(bp, 1000) # signal was sampled at 1000 Hz\nscatter!(peaks, bp[peaks], lab = \"Systolic pressure\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Cardio","category":"page"},{"location":"#Cardio.jl","page":"Home","title":"Cardio.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A simple toolbox that bundles functionality to work with cardiovascular time series. ","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"Utility.md\",\n    \"Detection.md\",\n    \"HRV.md\",\n    \"BRS.md\",\n]\nDepth = 1","category":"page"}]
}
