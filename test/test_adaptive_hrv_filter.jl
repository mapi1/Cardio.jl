#Test Values
test_values = Float64.(vec([852 845 846 851 846 846 840 825 823 821 836 854 854 850 822 805 802 804 804 823 815 799 802 793 792 807 821 823 809 809 823 840 858 847 823 799 788 788 789 799 805 805 817 853 873 896 900 886 862 840 845 854 858 861 830]))

# Ideal Case (No changes will be performed)
@testset "Ideal Case" begin
    fixed_signal, outliers, nonnormal = adaptive_hrv_filter(test_values)
    @test test_values == fixed_signal
    @test sum(outliers) == 0
    @test sum(nonnormal) == 0
end
fixed_signal, outliers, nonnormal = adaptive_hrv_filter(test_values)

# Changes but possibly physiological
changed_idxs = [20, 21]
changed_values = [630, 1340]
changed_test_values = copy(test_values)
changed_test_values[changed_idxs] = changed_values

@testset "Physiological Changes" begin
    fixed_signal, outliers, nonnormal = adaptive_hrv_filter(changed_test_values)
    notchanged = broadcast(!, nonnormal)
    @test sum(outliers) == 0
    @test changed_test_values[notchanged] == fixed_signal[notchanged]
end

# Changes to non physiological values
changed_idxs = [20, 21]
changed_values = [180, 2500]
changed_test_values = copy(test_values)
changed_test_values[changed_idxs] = changed_values

@testset "Non Physiological Changes" begin
    fixed_signal, outliers, nonnormal = adaptive_hrv_filter(changed_test_values)
    notchanged = broadcast(!, outliers)
    @test sum(nonnormal) == 0
    @test sum(outliers) == 2
    @test length(fixed_signal) + sum(outliers) == length(changed_test_values)
end

