TGv = v"4.0.3+TS"
# Generates template files to be used for spikesorting.

# * Changelog
#   Updated to TSHelperTS4.4.0
#   Updated def_template_gen_vars() name

# ! Errors
#   Sometimes ISIresult is missing a value

#   TODO List
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

pe = pyimport("matplotlib.patheffects")
# scisig = pyimport("scipy.signal")
# MultiCursor = matplotlib.widgets.MultiCursor
CheckButtons = matplotlib.widgets.CheckButtons

# Working directory and dependencies
cd(dirname(@__FILE__))
include("src/spike_sort_main.jl")
include("src/spike_sort_plotting.jl")
include("src/TSHelperTS4.4.0.jl")       #! Must remain after other dependencies

df_t = DataFrame(XLSX.readtable("savedFilesInfo2023v3.xlsx", "Descriptions"))
# df_t = df |>
#     @filter(_.Animal == animalID) |>
#     @filter(_.Position == pos) |>
#     @filter(_.Cell == cellDescriptor) |>
#DataFrame

for loadingoutput in eachrow(df_t)
    templateFilename = loadingoutput["Template"]
    animalID = loadingoutput["Animal"]
    stim_type = loadingoutput["Stimtype"]
    println("running next file: ", animalID, "  ", templateFilename)

    Random.seed!(42) #set the random seed for kmeans clustering, to allow for consistency

    # # * Parameters
    # # File Parameters
    # animalID = "SNT06"
    # pos = "1"
    # amp = "100"
    # block_amp = "0"
    # stim_type = "M1"
    # files_num = [x for x in range(1, length=5)]
    # channels = [5,14,29]
    k_override = 0

    # Spikesorting Parameters
    fs=30000

    if stim_type == "L1"
        pre_trial_len = round(Int, 1.0 * fs)
        post_trial_len = round(Int,5.0 * fs)
        WOI = [1000, 1500]
    elseif stim_type == "M3"
        pre_trial_len = round(Int, 0.5 * fs)
        post_trial_len = round(Int,1.5 * fs)
        WOI = [100, 600]
    elseif stim_type == "M2"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI = [0, 400]
    elseif stim_type == "M1"
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI = [0, 50]
    else
        pre_trial_len = round(Int, 0.4 * fs)
        post_trial_len = round(Int,0.8 * fs)
        WOI = [0, 100]
    end

    # WOI = [0, 100]
    onems = round(Int, fs/1000)
    # onems = isinteger(fs/1000) ? Int(fs/1000) : fs/1000
    pre_window = 7  # peak extraction
    post_window = 7 # peak extraction
    # pre_trial_len = round(Int, 0.4 * fs)
    # post_trial_len = round(Int,0.8 * fs)

    # PCA_ratio = 0.9
    # p_ratio = 0.95

    # MAX_CLUSTERS = 15
    # CLUSTERING_EPOCS = 1000
    # SHOW_CLUSTER_RESULT = true

    # Function Parameters
    # saveTemplate = false     # whether or not to save the template (ignored if you are loading from a template)
    # saveHistInspect = false     # whether or not to export raster data
    # autoMerge = true    # whether or not to merge all clusters
    # manualMerge = [[1, 2]]    # a vector of vectors indicating which clusters to merge - only used if autoMerge == false
    samples_on_x=false  # what units to use for x-axis of raster plots - true for samples or false for milliseconds
    showgraphs = false # toggle graph generation for faster debugging
    global ISIresult = [] #set up a global to store the ISI results


    ################################################################################



    templateFile = load("loadingOutput/$animalID/$templateFilename.jld2")
    signal_matrix, pk_bkgnds, file_lengths, stim_inds, cum_stim_inds = def_template_gen_vars(templateFile)

    # * template construction
    #' ## Peak Detection
    pk_inds, waveforms, _pk_heights = detect_spikes(signal_matrix,
                pk_bkgnds,
                stim_inds, 
                [pre_window, post_window],
                file_lengths = file_lengths, 
                window_of_interest = WOI .* onems)
                        
    #' ## Peak Clustering and Template Generation
    cluster_result, templates, PCA_Model, assign_order = cluster_spikes(waveforms,
                k_override = k_override,
                apply_pca=true,
                PCA_ratio = 0.99,
                PCA_scale = "off",
                max_clusters = 10,
                num_epocs = 100,
                SHOW_CLUSTER_RESULT = showgraphs)
    assigns = cluster_result.assignments
    println("Initial counts: ", counts(assigns), " (k-override = $k_override)")



    # * template application
    #' Get thresholds for clustering
    # use cluster_pratio = 1.0, then reduce until you remove the clearly unwanted waveforms. Conduct the following lines of code iteratively until satisfied
    thresholds = get_cluster_thresholds(waveforms, templates, assigns, cluster_pratio = 1.0, x_min = 10)
    # thresholds = get_cluster_thresholds(waveforms, templates, assigns, cluster_pratio = 0.95, x_min = 10)

    #' define data file to apply templates
    pk_inds_full, waveforms_full, _pk_heights_full = detect_spikes(signal_matrix, pk_bkgnds, stim_inds, [pre_window, post_window], file_lengths = file_lengths)
        
    #' Assigning spikes to templates
    assignment = assign_points_to_clusters(waveforms_full, templates, thresholds)

    #Plots of assigned waveforms
    if showgraphs
        pyplot_waveforms_separately(waveforms_full, assignment, centroids = templates, title = "Assigned from Templates")
    end

    # Combine all stim_inds
    pks_plot, ass_plot = assign_to_trials(cum_stim_inds, pk_inds_full, assignment, pre_trial_len, post_trial_len) 
    # pks_plot, ass_plot = assign_to_trials(cum_stim_inds, pk_inds2, assignment, 0, 300) 

    # pks for plotting, assignment for plotting
    #FIXME:     I've added pre/post trial lengths, but maybe these should go above somewhere, 
    #           e.g. the assignments above?
    #TODO:      Check I did this correctly: # cum_stim_inds = stim_inds2# Should replace for concat_stim_inds
    #           check we need the pre_trial_len and post_ are required

    # separate_events(pks_plot, ass_plot)
    ass_counts = counts(convert.(Int64, vcat(ass_plot...)))
    num_assigns = maximum(vcat(ass_plot...)) + 1
    numTrials = size(ass_plot, 1)
    pksLocs = Array{Any}(undef, maximum(vcat(ass_plot...))+1)
    pksAss = Array{Any}(undef, maximum(vcat(ass_plot...))+1)
    for i = 1:num_assigns
        pksLocs[i]=[]
        pksAss[i]=[]
        for j = 1:numTrials
            push!(pksLocs[i],pks_plot[j][findall(ass_plot[j] .== (i-1))])
            push!(pksAss[i],ass_plot[j][findall(ass_plot[j] .== (i-1))])
        end
        
        # ISI testing
        testSum = []
        for trial in pksLocs[i]
            push!(testSum, !isempty(trial))
        end
        if sum(testSum) > 0
            println("Assignment $i")
            result = isiTests(pksLocs[i], onems)
            println(result, "\n")
            push!(ISIresult, result)
        else println("Assignment $i is empty\n")
            push!(ISIresult, "EMPTY")
        end
    end

    Mv = v"4.1.1+FPA"
    # For modifying templates after their initial generation. Separated from templateGen for convenience.
    # Run directly after templateGen (i.e. without killing REPL). Note that part of this code is uncommented for easier selection - refer to templateGen for commenting.

    # * Changelog
    #   Added fully automated loop to retemplate assignments that fail the ISI test
    #   Added ISIresult to the saved output

    # ! Errors

    #   TODO List
    # figure out how to iterate assigns, templates, templatecounternew - it doesn't like replacing them

    # ? Suggestions & Queries


    global retemplating = true #to initiate the function
    global k_override = 0 # ! this will override every looped retemplate, so should be deprecated

    ## TODO check here for templating to make sure we don't do it before we need to

    global templateinit = zeros(length(unique(assigns))) #make a counter to store how many retemplates have occured
    global retemplate = 1
    global templatecounter = []
    global templatecounternew = []

    while retemplating

        # redefine variables from templategen so that they have global scope and can be parsed in the while loop
        # global assigns = assigns
        # global templates = templates
        # global thresholds = thresholds
        # global pk_inds_full = pk_inds_full
        # global waveforms_full = waveforms_full
        # global _pk_heights_full = _pk_heights_full
        # global assignment = assignment
        # global pks_plot = pks_plot
        # global ass_plot = ass_plot
        # global ass_counts = ass_counts
        # global num_assigns = num_assigns
        # global numTrials = numTrials
        # global pksLocs = pksLocs
        # global pksAss = pksAss 
        # global templateinit = templateinit
        # global retemplate = retemplate
        # global templatecounter = templatecounter
        # global templatecounternew = templatecounternew

        if templatecounter == [] 
            templatecounter = templateinit
        else 
            templatecounter = templatecounternew
        end

        if retemplate > length(templatecounter)
            println("All retemplating completed!")
            break
        end

        if templatecounter[retemplate] > 3 || ISIresult[retemplate + 1] != "FAILED"
            retemplate += 1
            println("Templating limit reached for this assignment, moving to next assignment")
            continue

        else
            showgraph = false

            #retemplate -= 1     #! DO NOT REMOVE - this transforms retemplate into the assignment number before unassigned is added
            waveformsRetemplate = fndAssign(waveforms, retemplate, assigns)



            # Retemplate single assignment
            cluster_resultNew, templatesNew, PCA_ModelNew, assign_orderNew = cluster_spikes(waveformsRetemplate,
                        k_override=k_override,
                        apply_pca=true,
                        PCA_ratio=0.99,
                        PCA_scale="off",
                        max_clusters=10,
                        num_epocs=100,
                        SHOW_CLUSTER_RESULT=showgraph)
            assignsNew = cluster_resultNew.assignments
            println("Retemplated counts: ", counts(assignsNew), " (k-override = $k_override)")



            """
            # reorder templates to correct assignments if required and rerun figures above to check assignments: 
            templates = [templates[:,3] templates[:,1] templates[:,2] templates[:,4] templates[:,5]]
            counts(assigns) # see cluster numbers and check order. They should be descending if your templates were correctly ordered
            # or:
            # templates = templates[:,assign_order] # not working atm, use approach on line above

            # Other tools to help you ensure good centroids:
                # use the following tools only if you cannot get clean clusters from above. Sometimes you need to remove a dominant cluster then recluster to get the correct centoids

                # to keep assignments:
                waveforms1=fndAssign(waveforms, 1, assigns) # this only keeps waveforms from assignment 1

                # to remove assignments: 
                waveforms1=rmAssign(waveforms, 1, assigns) # this keeps all waveforms except assignment 1

                #use this for merging clusters. It might be better to do it manually as above with the fndAssign() and rmAssign() and curate accordingly, that gives you more control over the final result. The merging in theory can be avoided entirely if the cluster result works ok. 
                merges = autoMerge ? [[n] for n in 1:maximum(assigns)] : manualMerge #! use this as last resort, try to get your clustering solutions correct from the beginning
                #templates = keep_and_merge_templates(templs, merges)
                merges = [[1, 2], [3]]
                templates, assigns = keep_and_merge_templates(templates,merges,assigns)
                
            """

            assignsOriginal = assigns
            templatesOriginal = templates
            global ISIresult = []

            assigns, templates, templatecounternew = combine_templates(retemplate, assignsOriginal, templatesOriginal, assignsNew, templatesNew, templatecounter, showgraph)
            thresholds = get_cluster_thresholds(waveforms, templates, assigns, cluster_pratio = 1.0, x_min = 10)
            pk_inds_full, waveforms_full, _pk_heights_full = detect_spikes(signal_matrix, pk_bkgnds, stim_inds, [pre_window, post_window], file_lengths = file_lengths)
            assignment_temp = assign_points_to_clusters(waveforms, templates, thresholds)
            assignment = assign_points_to_clusters(waveforms_full, templates, thresholds)
            if showgraph
                pyplot_waveforms_separately(waveforms_full, assignment, centroids = templates, title = "Assigned from Templates")
            end
            pks_plot, ass_plot = assign_to_trials(cum_stim_inds, pk_inds_full, assignment, pre_trial_len, post_trial_len)
            ass_counts = counts(convert.(Int64, vcat(ass_plot...)))
            num_assigns = maximum(vcat(ass_plot...)) + 1
            numTrials = size(ass_plot, 1)
            pksLocs = Array{Any}(undef, maximum(vcat(ass_plot...))+1)
            pksAss = Array{Any}(undef, maximum(vcat(ass_plot...))+1)
            for i = 1:num_assigns
                pksLocs[i]=[]
                pksAss[i]=[]
                for j = 1:numTrials
                    push!(pksLocs[i],pks_plot[j][findall(ass_plot[j] .== (i-1))])
                    push!(pksAss[i],ass_plot[j][findall(ass_plot[j] .== (i-1))])
                end
                testSum = []
                for trial in pksLocs[i]
                    push!(testSum, !isempty(trial))
                end
                if sum(testSum) > 0
                    println("Assignment $i")
                    result = isiTests(pksLocs[i], onems)
                    println(result, "\n")
                    push!(ISIresult, result)
                else println("Assignment $i is empty\n")
                    push!(ISIresult, "EMPTY")
            end
        end
        println("Retemplating completed, checking next assignment...")
        end
    end


    #* Plotting
    # pyplot_waveforms_separately(waveforms_full, assignment, centroids = templates, title = "Assigned from Templates")
    # for waveNo in 1:length(templates[1,:])
    #     temporaryWave = fndAssign(waveforms, waveNo, assigns)
    #     figure()
    #     PyPlot.plot(temporaryWave, color=COLOURS[waveNo])
    #     PyPlot.plot(templates[:, waveNo], color="black")
    # end
    # pks_plot_zeroed = Vector()
    # for pk_trial in pks_plot
    #     pk_trial_zeroed = samples_on_x ? pk_trial .- pre_trial_len : (pk_trial .- pre_trial_len) ./ onems # subtracts the pre-trial length and converts to milliseconds if specified
    #     push!(pks_plot_zeroed, pk_trial_zeroed)
    # end
    # pyplot_rasters(pks_plot_zeroed, ass_plot, samples_on_x, pre_trial_len, post_trial_len, title="Raster of assignments")
    # for i = 1:num_assigns
    #     pksLocs_zeroed = Vector()
    #     for pk_trial in pksLocs[i]
    #         pk_trial_zeroed = samples_on_x ? pk_trial .- pre_trial_len : (pk_trial .- pre_trial_len) ./ onems # subtracts the pre-trial length and converts to milliseconds if specified
    #         push!(pksLocs_zeroed, pk_trial_zeroed)
    #     end
    #     testSum = []
    #     for trial in pksLocs[i]
    #         push!(testSum, !isempty(trial))
    #     end
    #     if sum(testSum) > 0
    #         pyplot_rasters(pksLocs_zeroed, pksAss[i], samples_on_x, pre_trial_len, post_trial_len, title="Raster of assignments = $i")
    #     end
    # end
    #* Save output
    #mkdir"spikesortingTemplate/$animalID"
    mkpath("spikesortingTemplate/$animalID")
    @save "spikesortingTemplate/$animalID/$templateFilename.jld2" waveforms templates assigns thresholds ISIresult
    println("Template generation complete, moving to next file..")
end
println("Batch templating complete!")