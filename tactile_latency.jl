TLv = v"4.0.0+TS"
# Exports traces and spike data as CSVs for plotting.

# * Changelog

# ! Errors

#   TODO List

# ? Suggestions & Queries


# Packages
using CSV
using DataFrames
using FileIO
using JLD2
using Query
using XLSX

# Working directory and dependencies
cd(dirname(@__FILE__))


bins = 0:49
units_list_df = CSV.read("analysisResults/tactile_evoked_units.csv", DataFrame)
counts_df = DataFrame((NeuronID=[], Treatment=[], PeakLatency=[]))

for unit in eachrow(units_list_df)
    animalID = unit["animal_ID"]
    pos = unit["position"]
    cell_desc = unit["cell_desc"]
    condition = unit["stim_type"]
    comparison_files = unit["recording_files_num"]
    template_files = unit["template_files_num"]
    channels = unit["channels"]
    comparisonFilename = "$pos 100 0 M1 $comparison_files $channels"
    templateFilename = "$pos 100 0 M1 $template_files $channels"

    filename = string(animalID*"-"*string(pos)*"-"*cell_desc*" "*condition*"-"*comparisonFilename*"-"*templateFilename)
    ssO_file = load("spikesortingOutput/$filename.jld2")


    assignment = unit["assignment"]
    pks_plot = ssO_file["pks_plot"]
    ass_plot = ssO_file["ass_plot"]
    trials_count = ssO_file["trials_count"]
    pre_trial_len = ssO_file["pre_trial_len"]

    pksLocs = []
    for j = 1:trials_count
        push!(pksLocs, pks_plot[j][findall(ass_plot[j] .== (assignment-1))])
    end


    raster_df = DataFrame((Time=[], Trial=[]))
    for (trial_no, pk_trial) in enumerate(pksLocs)
        for pk in pk_trial
            pk_zeroed = (pk - pre_trial_len) / 30 # subtracts the pre-trial length and converts to milliseconds
            if 0 <= pk_zeroed <= 50
                push!(raster_df, [pk_zeroed, trial_no])
            end
        end
    end

    (bin_count, latency) = (0, 0)
    for bin_start in bins
        bin_end = bin_start + 1

        bin_count_new = count(i->(bin_start<i<=bin_end), raster_df[!, "Time"])

        if bin_count_new > bin_count
            (bin_count, latency) = (bin_count_new, bin_start)
        end
    end

    
    neuronID = unit["neuron_ID"]
    treatment = unit["treatment"]
    push!(counts_df, (neuronID, treatment, latency))

    println(neuronID)
end


CSV.write("analysisResults/tactile_peak_latencies.csv", counts_df)
