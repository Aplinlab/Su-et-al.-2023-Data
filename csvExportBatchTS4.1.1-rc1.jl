CEbv = v"4.1.1-rc1+TS"
# Exports spikesortingOutput to a CSV after counting spikes.

# * Changelog
#   Added std to output
#   Commented out redundant(?) code defining result
#!      Needs to be tested

# ! Errors
#   Empty assignment workaround may occasionally produce an error (needs more documentation)

#   TODO List
#   Individual trial counts
#   Better time information (pull time from filename according to file number)
#   Clean up for upload

# ? Suggestions & Queries



# Packages
using DSP
using Glob
using LinearAlgebra
using MultivariateStats
using Printf
using PyPlot; pygui(true) # for Display

using Statistics #std
using StatsBase #mode

using CSV
using JLD2
using FileIO
using PyCall
using DataFrames
using XLSX
using Query
using Dates
# scisig = pyimport("scipy.signal")
# MultiCursor = matplotlib.widgets.MultiCursor
CheckButtons = matplotlib.widgets.CheckButtons

# Working directory and dependencies
cd(dirname(@__FILE__))
include("src/spike_sort_main.jl")
include("src/spike_sort_plotting.jl")
include("src/TSHelperTS4.4.0.jl")       #! Must remain after other dependencies



# * Parameters
# Function Parameters
ignore_unassigned = true
samples_on_x=false  # what units to use for x-axis of raster plots - true for samples or false for milliseconds



################################################################################

df = DataFrame(XLSX.readtable("savedFilesInfo2023v3.1.xlsx", "lOutput and sTemplate"))

outputTable = DataFrame(
    neuron_ID = [], 

    animal_ID = [], 
    treatment = [], 
    sex = [], 
    position = [], 
    channels = [], 
    cell_desc = [], 
    assignment = [], 
    cell_type = [], 
    stim_type = [], 
    stim_code = [], 
    stim_amp = [], 
    idc_amp = [], 
    recovery_time = [], 
    exp_phase = [], 
    is_good_assign = [], 

    WOI1 = [], WOI2 = [], 
    spikes_count_pre = [], spikes_count_post = [], trials_count = [], 
    spikes_std_pre = [], spikes_std_post = [], 
    recording_time = [], template_time = [], delta_t = [], 

    spikes_pre = [], spikes_post = [], 
    max_bin_time = [], max_bin_count = [], max_bin_time_template = [], max_bin_count_template = [], bins = [], 

    recording_filenames = [], recording_files_num = [], 
    template_filenames = [], template_files_num = [], 
)

date = Dates.today()

for comparison in eachrow(df)
    templateFilename = comparison["Template"]
    comparisonFilename = comparison["Filename"]
    animalID = comparison["Animal"]
    pos = parse(Int, comparison["Position"])
    cell_desc = comparison["Cell"]
    condition = comparison["Condition"]
    recovery_time = parse(Float64, comparison["Recovery Time"])

    println("$animalID $pos $cell_desc $condition ($comparisonFilename)")

    description = string(animalID*"-"*string(pos)*"-"*cell_desc*" "*condition*"-"*comparisonFilename*"-"*templateFilename)
    spikesortedFile = load("spikesortingOutput/$description.jld2")

    animalID_file, treatment, sex, pos_file, channels, cell_desc_file, cell_type, stim_type, stim_code, stim_amp, idc_amp, recovery_time_file, exp_phase,
        WOI1, WOI2, trials_count, recording_time, template_time, delta_t, recording_filenames, recording_files_num, template_filenames, template_files_num,
        pks_plot, ass_plot, pre_trial_len, post_trial_len, onems, pre_window, post_window, num_assigns, good_assigns = def_csv_vars(spikesortedFile)
    
    # Error checking
    if animalID != animalID_file
        error("Animal IDs do not match")
    elseif pos != pos_file
        error("Positions do not match")
    elseif cell_desc != cell_desc_file
        error("Cell descriptions do not match")
    elseif recovery_time != recovery_time_file
        error("Recovery times do not match")
    end

    # Overwriting sex for SPC18
    if animalID == "SPC18"
        sex = "F"
    end

    # Overwriting WOI for pinch cells
    if stim_code == "M3"
        WOI2 = [50, 400]
        WOI1 = [WOI2[1] - WOI2[2], 0]
    end

    # Define bin sizes
    if stim_code == "L1"
        bins = WOI2[1]:50:WOI2[2]
    elseif stim_code == "M3"
        bins = WOI2[1]:10:WOI2[2]
    elseif stim_code == "M2"
        bins = WOI2[1]:1:WOI2[2]
    elseif stim_code == "M1"
        bins = WOI2[1]:1:WOI2[2]
    else
        bins = WOI2[1]:1:WOI2[2]
    end

    ################################################################################

    # separate_events(pks_plot, ass_plot)
    ass_counts = counts(convert.(Int64, vcat(ass_plot...)))
    pksLocs = Array{Any}(undef, num_assigns)
    pksAss = Array{Any}(undef, num_assigns)
    first_assignment = ignore_unassigned ? 2 : 1
    for i = first_assignment:num_assigns
        if i > length(ass_counts)
            push!(ass_counts, 0)
        end
        pksLocs[i]=[]
        pksAss[i]=[]
        for j = 1:trials_count
        push!(pksLocs[i],pks_plot[j][findall(ass_plot[j] .== (i-1))])
        push!(pksAss[i],ass_plot[j][findall(ass_plot[j] .== (i-1))])
        end

        testSum = []
        for trial in pksLocs[i]
            push!(testSum, !isempty(trial))
        end
        
        pksLocs_zeroed = Vector()
        for pk_trial in pksLocs[i]
            pk_trial_zeroed = samples_on_x ? pk_trial .- pre_trial_len : (pk_trial .- pre_trial_len) ./ onems # subtracts the pre-trial length and converts to milliseconds if specified
            push!(pksLocs_zeroed, pk_trial_zeroed)
        end
        
        # if sum(testSum) > 0
        #     # pyplot_rasters(pksLocs_zeroed, pksAss[i], samples_on_x, pre_trial_len, post_trial_len, title="Raster of assignments = $i")
        #     # println("Assignment $i")
        #     result = isiTests(pksLocs[i], onems, reject=false)
        #     println(result, "\n")
        # else
        #     result = "PASSED"
        #     println("Assignment $i is empty\n")
        # end

        neuronID = string(animalID*"_"*string(pos)*"_"*string(channels')*"-"*string(i)*"_"*cell_desc)
        spikes_pre = []
        spikes_post = []
        for trial in pksLocs_zeroed
            for spike in trial
                if WOI1[1] < spike <= WOI1[2]
                    push!(spikes_pre, spike)
                elseif WOI2[1] < spike <= WOI2[2]
                    push!(spikes_post, spike)
                end
            end
        end
        spikes_count_pre_pertrial = [count(a->(WOI1[1]<a<=WOI1[2]), x) for x in pksLocs_zeroed]
        spikes_count_post_pertrial = [count(a->(WOI2[1]<a<=WOI2[2]), x) for x in pksLocs_zeroed]
        spikes_count_pre = sum(spikes_count_pre_pertrial)
        spikes_count_post = sum(spikes_count_post_pertrial)
        if length(spikes_pre) != spikes_count_pre
            error("Pre-stim counts don't match")
        elseif length(spikes_post) != spikes_count_post
            error("Post-stim counts don't match")
        end
        spikes_std_pre = Statistics.std(spikes_count_pre_pertrial)
        spikes_std_post = Statistics.std(spikes_count_post_pertrial)

        max_bin = [0, 0]
        for bin_start in bins
            bin_end = bin_start + bins.step
            bin_count = sum(count(a->(bin_start<a<=bin_end), x) for x in pksLocs_zeroed)
            if bin_count > max_bin[2]
                max_bin = [bin_start, bin_count]
            end
        end

        result = isiTests(pksLocs[i], onems, reject=false)
        if i > length(good_assigns)
            is_good_assign = "There are too few items in ISI results (unknown bug in templateGen). Investigate $animalID $pos $cell_desc further."
        else is_good_assign = result == "PASSED" && good_assigns[i] == true
        end
        push!(outputTable, (
            neuronID, 

            animalID, 
            treatment, 
            sex, 
            pos, 
            channels, 
            cell_desc, 
            i, 
            cell_type, 
            stim_type, 
            stim_code, 
            stim_amp, 
            idc_amp, 
            recovery_time, 
            exp_phase, 
            is_good_assign, 

            WOI1, WOI2, 
            spikes_count_pre, spikes_count_post, trials_count, 
            spikes_std_pre, spikes_std_post, 
            recording_time, template_time, delta_t, 

            spikes_pre, spikes_post, 
            max_bin[1], max_bin[2], 0, 0, bins, 

            recording_filenames, recording_files_num, 
            template_filenames, template_files_num, 
        ))
    end
end

# for (row_no, unit) in enumerate(eachrow(outputTable))
#     template_row = outputTable |>
#         @filter(_.recording_filenames == unit["template_filenames"]) |>
#         @filter(_.neuron_ID == unit["neuron_ID"]) |>
#     DataFrame

#     if size(template_row)[1] == 1
#         unit["max_bin_time_template"] = template_row[1, "max_bin_time"]
#         unit["max_bin_count_template"] = template_row[1, "max_bin_count"]
#     elseif size(template_row)[1] > 1
#         error("Too many matching template rows (row $row_no)")
#     else error("No matching template rows (row $row_no)")
#     end

#     println(row_no)
# end

CSV.write("analysisResults/fullResults $date.csv", outputTable)
println("CSV Exporting Complete")