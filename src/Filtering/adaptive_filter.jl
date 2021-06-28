"""
adaptive_hrv_filter(signal::Vector{<:Real}; remove_outliers::Bool = true, replace_nonnormal::Bool = true, adaptive_controlling_coef::Real = 0.05, proportionality_limit::Real = 10/100, outlier_minz_factor::Real = 3, max_excess_hrv::Real = 20, physiological_values::Tuple{Real, Real} = (200, 2000))

An addaptive filter for HRV data
     Based on: Wessel, N., Voss, A., Malberg, H., Ziehmann, Ch., 
        Voss, H. U., Schirdewan, A., Meyerfeldt, U.,Kurths, J.: 
        Nonlinear analysis of complex phenomena in cardiological data, 
        Herzschr. Elektrophys., 11(3), 2000, 159-173, doi:10.1007/s003990070035.

Filters out unphysiological beats in the RR series with the ability to replace them.

# Args:

* 'signal::Vector': HRV in ms

# Keyword Args
* 'remove_outliers::Bool = true': option if non physiological outliers shall be removed
* 'replace_nonnormal::Bool = true': if true, non normal HRV values are replaced
* 'adaptive_controlling_coef::Real = 0.05': ???
* 'proportionality_limit::Real = 10/100': ???
* 'outlier_minz_factor::Real = 3': ???
* 'max_excess_hrv::Real = 20': ???
* 'physiological_values::Tuple{Real, Real} = (200, 2000)': Definition of the physiological values, scheme: (min, max)

# Return
Return value depends on keyword args

```julia
julia> adaptive_hrv_filter(signal)
```
"""
function adaptive_hrv_filter(signal::Vector{<:Real}; remove_outliers::Bool = true, replace_nonnormal::Bool = true, adaptive_controlling_coef::Real = 0.05, proportionality_limit::Real = 10/100, outlier_minz_factor::Real = 3, max_excess_hrv::Real = 20, physiological_values::Tuple{Real, Real} = (200, 2000))
    # some input checking
    0 < adaptive_controlling_coef < 1 || throw(DomainError("adaptive_controlling_coef has to be in range 0:1."))
    physiological_values[1] < physiological_values[2] || throw(DomainError("physiological_values is a tuple with minumum first and maximum second, is: $physiological_values."))
    # TODO

    # Remove values that are out of range of Physiological Heart beats
    if remove_outliers
        ok = (physiological_values[1] .<= signal .<= physiological_values[2])
        removed_nonphysiological_outliers = broadcast(!, ok)
        signal = signal[ok]
        @assert length(signal) >= 1 "Signal contains only outliers, recheck parameters"
    end
    # Calculate helper matrixes for the Binominal Filter
    bin_coeff = [1; 6; 15; 20; 15; 6; 1]
    coeff_sum = sum(bin_coeff)
    # normalize filter coefs
    bin_coeff =  bin_coeff .* (1/coeff_sum)
    coeff_count = length(bin_coeff)
    @assert (mod(coeff_count, 2) == 1) "Nono"
    
    # signalLen = length(signal)
    # filter_values = repeat(bin_coeff', signalLen)'
    # value_index_column = (1:coeff_count) .- Int(ceil(coeff_count/2))
    # hrv_ind_offset_val = repeat(value_index_column', signalLen)
    # hrv_ind_base_val = repeat(collect(1:signalLen)', coeff_count)'
    
    # #The first and last elements will be out of range work around by just (re-)using first and last elements...
    # hrv_idx_mat = hrv_ind_offset_val .+ hrv_ind_base_val
    # replace!(x -> x <=0 ? 1 : x, hrv_idx_mat)
    # replace!(x -> x > signalLen ? signalLen : x, hrv_idx_mat)

    # Calculate first filtered Signal (Through Binominal Filter)
    # filtered_signal = vec(sum(signal[hrv_idx_mat] .* filter_values', dims = 2)) 
    filtered_signal = filtfilt(bin_coeff, signal)

    # Apply "Adaptive Percent Filter" => fixed_hrv_timeseries
    adaptive_mean, adaptive_sigma = adaptive_moments(filtered_signal, coeff_count, adaptive_controlling_coef)
    adaptive_sigma_mean = mean(adaptive_sigma)

    last_good_value = signal[1] #Maybe filtered_signal(1)
    last_good_range = proportionality_limit * last_good_value + outlier_minz_factor * adaptive_sigma_mean

    hrv_diff = diff(signal)
    normal_values = repeat([true], length(filtered_signal)) 
    for i in 2: length(signal) 
        current_diff = abs(hrv_diff[i - 1])
        current_max_range = proportionality_limit * signal[i-1] + outlier_minz_factor * adaptive_sigma_mean
        current_value_is_normal= current_diff <= current_max_range || current_diff <= last_good_range
        if current_value_is_normal
            last_good_value = signal[i]
            last_good_range = proportionality_limit * last_good_value + outlier_minz_factor * adaptive_sigma_mean
        end
        normal_values[i] = current_value_is_normal 
    end
    # Just change reference, no copy
    fixed_signal = signal
    nonnormal_values = broadcast(!, normal_values)
    fixed_signal[nonnormal_values] = adaptive_mean[nonnormal_values] + (randn(sum(nonnormal_values)) .- 0.5) .* adaptive_sigma[nonnormal_values]
    
    # Calculate second filtered Signal (Through Binominal Filter)
    #filtered_fixed_signal =  vec(sum(fixed_signal[hrv_idx_mat] .* filter_values', dims = 2))
    filtered_fixed_signal = filtfilt(bin_coeff, fixed_signal) 

    # Apply "Adaptive Controlling Procedure" => fixed_fixed_signal
    adaptive_fixed_mean, adaptive_fixed_sigma = adaptive_moments(filtered_fixed_signal, coeff_count, adaptive_controlling_coef)

    normal_fixed_hrv_values = abs.(fixed_signal - adaptive_fixed_mean) .<= (outlier_minz_factor .* adaptive_fixed_sigma .+ max_excess_hrv)

    fixed_fixed_signal = fixed_signal
    nonnormal_fixed_values = broadcast(!, normal_fixed_hrv_values)
    fixed_fixed_signal[nonnormal_fixed_values] = adaptive_fixed_mean[nonnormal_fixed_values] + (randn(sum(nonnormal_fixed_values)) .- 0.5) .* adaptive_fixed_sigma[nonnormal_fixed_values]

    # Returns
    if replace_nonnormal
        if remove_outliers
            return fixed_fixed_signal, removed_nonphysiological_outliers, nonnormal_values 
        else
            return fixed_fixed_signal, nonnormal_values 
        end
    else
        if remove_outliers
            return removed_nonphysiological_outliers, nonnormal_values 
        else
            return nonnormal_values 
        end
    end
end

# Helper to calculate adaptive moments (mean, std)
function adaptive_moments(signal::Vector{<:Real}, init_length::Int, adaptive_controlling_coef::Real)
    adaptive_mean = Float64[]
    adaptive_variance = Float64[]
    push!(adaptive_mean, mean(signal[1:init_length]))
    push!(adaptive_variance, var(signal[1:init_length]))
    for i in 2:length(signal)
        push!(adaptive_mean, (adaptive_mean[end] - adaptive_controlling_coef * (adaptive_mean[end] - signal[i-1])))
        last_variance_item = (adaptive_mean[end-1] - signal[i-1]) ^ 2
        push!(adaptive_variance, adaptive_variance[end] - adaptive_controlling_coef * (adaptive_variance[end] - last_variance_item))
    end  
   return adaptive_mean, sqrt.(adaptive_variance)
end
