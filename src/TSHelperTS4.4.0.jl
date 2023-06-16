TSHv = v"4.4.0+TS"
# Helper functions and constants for most of the spikesorting pipeline.

# * Changelog
#   Added def_csv_vars()
#   Shortened variable in function names to var

# ! Errors

#   TODO List

# ? Suggestions & Queries



using PyCall

COLOURS = ["#e69f00ff", "#56b4e9ff", "#009e73ff", "#0072b2ff", "#d55e00ff", "#cc79a7ff", "#f0e442ff", "gray", "olive", "orange", "#8c564b", "#e377c2", "#e69f00ff", "#56b4e9ff", "#009e73ff", "#0072b2ff", "#d55e00ff", "#cc79a7ff", "#f0e442ff", "gray", "olive", "orange", "#8c564b", "#e377c2"]
colours = ["#e69f00ff", "#56b4e9ff", "#009e73ff", "#0072b2ff", "#d55e00ff", "#cc79a7ff", "#f0e442ff", "gray", "olive", "orange", "#8c564b", "#e377c2", "#e69f00ff", "#56b4e9ff", "#009e73ff", "#0072b2ff", "#d55e00ff", "#cc79a7ff", "#f0e442ff", "gray", "olive", "orange", "#8c564b", "#e377c2"]

Random.seed!(42) #set the random seed for kmeans clustering, to allow for consistency

function pyplot_add_centroids(centroids)
    D,K = size(centroids)
    # css = Plots.palette(:darkrainbow,max(K,2)).colors.colors
    # colours = [(p.r, p.g, p.b) for p in css]
    global colours
    for k = 1:K
        PyPlot.plot(centroids[:,k], color = colours[k], linewidth=2, path_effects = [pe.Stroke(linewidth = 5, foreground = "black"), pe.Normal()])
    end
end

function pyplot_waveforms_separately(waveform, assignments, rows=3, cols=1; alpha=0.5, centroids=nothing, title="")
    global colours #assume enough colours
    D,N = size(waveform)
    K = maximum(assignments)
    rows = K
    num_plots = rows * cols
    f,ax = PyPlot.subplots(rows, cols, sharex=true, sharey=false)
    if K == 1 # Changed @17/02/2021 to allow single template plotting
        ax = [ax]
    end
    f.suptitle(title)
    for p = 1:min(K,num_plots)
        ax[p].plot(waveform[:,p .== assignments], color=colours[p], alpha=alpha)
        # if isnothing(centroids) #! removed if statement
            ax[p].plot(centroids[:,p], color="black")
        # end
    end
end

function pyplot_rasters(trials_of_pks, trials_of_pk_assignments, samples, pre_trial_len, post_trial_len; axes=nothing, tophigh=true, title="", stim_index=0, fs=nothing, visible_option=false)
    global colours
    PyPlot.ion()
    num_assigns = 1+maximum(vcat(trials_of_pk_assignments...)) #includes and assume there always is a zero assignment
    rasters = [[] for n in range(0,length=num_assigns)] # num_assigns length of rasters
    if isnothing(axes)
        f, axes = PyPlot.subplots(1,1)
    end
    num_trials = length(trials_of_pks)
    # trial_id is the trial number, pk_trial is the peak indices in the trial
    for (trial_id, pk_trial) in enumerate(trials_of_pks)
        local num_pks = length(pk_trial)
        pk_y  = tophigh ? trial_id : (num_trials-trial_id)
        for (pk_id, pk_x) in enumerate(pk_trial)
            pk_ass = trials_of_pk_assignments[trial_id][pk_id]
            c = pk_ass > 0 ? colours[pk_ass] : "lightgray"
            raster_dot = axes.plot(pk_x, pk_y, ".", color=c)[1]
            push!(rasters[pk_ass+1], raster_dot)
        end
    end
    axes.set_title(title)

    if visible_option
        labels = [string(n) for n in range(0,length=num_assigns)]
        visible = [true for n in range(0,length=num_assigns)]
        checkButton_region = PyPlot.axes([0.03, 0.4, 0.15, 0.15]) #axes copied from @Jason's code spikesorting_template_example_jp1.jl:253
        chxbox = CheckButtons(checkButton_region, labels, visible)

        function set_visible(label)
            # index = indexin([label],labels)[1]
            index = findall(x -> x==label, labels)[1] # Probably more neater indexing than above @Jason what do you think
            for raster_dot in rasters[index]
                raster_dot.set_visible(! raster_dot.get_visible())
            end
            PyPlot.draw() # NOT SURE IF THIS IS NEEDED
        end
        chxbox.on_clicked(set_visible)
    end

    PyPlot.vlines(0, 0, length(trials_of_pks)+1, colors="black" )
    PyPlot.ylabel("Trial no.")
    margin = (pre_trial_len+post_trial_len)/25
    if samples
        PyPlot.xlabel("Samples (collected at $onems kHz)")
        PyPlot.xlim(-pre_trial_len-margin,post_trial_len+margin)
    else
        PyPlot.xlabel("Time (ms)")
        PyPlot.xlim((-pre_trial_len-margin)/onems,(post_trial_len+margin)/onems)
    end
end

function def_template_gen_vars(wave_data::Dict{String, Any})
    signal_matrix = wave_data["signal_matrix"]
    pk_bkgnds = wave_data["pk_bkgnds"]
    file_lengths = wave_data["file_lengths"]
    stim_inds = wave_data["stim_inds"]
    cum_stim_inds = wave_data["cum_stim_inds"]

    return signal_matrix, pk_bkgnds, file_lengths, stim_inds, cum_stim_inds
    println("Variables extracted")
end

function def_spikesorting_vars(wave_data::Dict{String, Any}, template_data::Dict{String, Any}, template_file::Dict{String, Any})
    wave_signal_matrix = wave_data["signal_matrix"]
    wave_pk_bkgnds = wave_data["pk_bkgnds"]
    wave_file_lengths = wave_data["file_lengths"]
    wave_stim_inds = wave_data["stim_inds"]
    wave_cum_stim_inds = wave_data["cum_stim_inds"]

    template_waveforms = template_file["waveforms"]
    template_templates = template_file["templates"]
    template_assigns = template_file["assigns"]
    template_thresholds = template_file["thresholds"]
    template_isi_results = template_file["ISIresult"]

    stim_amp = parse(Float64, wave_data["amp"])
    idc_amp = parse(Int, wave_data["block_amp"])
    stim_code = wave_data["stim_type"]
    channels = wave_data["channels"]
    filenames = wave_data["filenames"]
    files_num = wave_data["files_num"]
    template_filenames = template_data["filenames"]
    template_files_num = template_data["files_num"]
    good_assigns = [x=="PASSED" for x in template_isi_results]

    return wave_signal_matrix, wave_pk_bkgnds, wave_file_lengths, wave_stim_inds, wave_cum_stim_inds, template_waveforms, template_templates, template_assigns, template_thresholds, stim_amp, idc_amp, stim_code, channels, filenames, files_num, template_filenames, template_files_num, good_assigns
    println("Variables extracted")
end

function def_csv_vars(spikesorted_data::Dict{String, Any})
    animalID = spikesorted_data["animalID"]
    treatment = spikesorted_data["treatment"]
    sex = spikesorted_data["sex"]
    pos = spikesorted_data["pos"]
    channels = spikesorted_data["channels"]
    cell_desc = spikesorted_data["cell_desc"]
    cell_type = spikesorted_data["cell_type"]
    stim_type = spikesorted_data["stim_type"]
    stim_code = spikesorted_data["stim_code"]
    stim_amp = spikesorted_data["stim_amp"]
    idc_amp = spikesorted_data["idc_amp"]
    recovery_time = spikesorted_data["recovery_time"]
    exp_phase = spikesorted_data["exp_phase"]
    WOI1 = spikesorted_data["WOI1"]
    WOI2 = spikesorted_data["WOI2"]
    trials_count = spikesorted_data["trials_count"]
    recording_time = spikesorted_data["recording_time"]
    template_time = spikesorted_data["template_time"]
    delta_t = spikesorted_data["delta_t"]
    recording_filenames = spikesorted_data["recording_filenames"]
    recording_files_num = spikesorted_data["recording_files_num"]
    template_filenames = spikesorted_data["template_filenames"]
    template_files_num = spikesorted_data["template_files_num"]
    pks_plot = spikesorted_data["pks_plot"]
    ass_plot = spikesorted_data["ass_plot"]
    pre_trial_len = spikesorted_data["pre_trial_len"]
    post_trial_len = spikesorted_data["post_trial_len"]
    onems = spikesorted_data["onems"]
    pre_window = spikesorted_data["pre_window"]
    post_window = spikesorted_data["post_window"]
    num_assigns = spikesorted_data["num_assigns"]
    good_assigns = spikesorted_data["good_assigns"]

    return animalID, treatment, sex, pos, channels, cell_desc, cell_type, stim_type, stim_code, stim_amp, idc_amp, recovery_time, exp_phase, WOI1, WOI2, trials_count, recording_time, template_time, delta_t, recording_filenames, recording_files_num, template_filenames, template_files_num, pks_plot, ass_plot, pre_trial_len, post_trial_len, onems, pre_window, post_window, num_assigns, good_assigns
    println("Variables extracted")
end

function detect_spikes(
    signal_matrix::AbstractMatrix, 
    pk_bkgnds::AbstractVector, 
    stim_inds::AbstractVector, 
    extract_window = [10,10]; 
    file_lengths::AbstractVector = [size(signal_matrix, 2)], 
    window_of_interest = nothing,
    intra_chan_window=6,
    inter_chan_window=6
    )
    
    num_channels = size(signal_matrix,1)
    sig_len = size(signal_matrix, 2)
    # onems = round(Int, fs/1000) # WOI provided in units of sample space

    largest_peaks = zeros(sig_len+1)
    cum_len = 1
    for (fn,fl) in enumerate(file_lengths)
        sis = stim_inds[fn]
        for ch_ind = 1:num_channels
            for si in sis
                si += cum_len
                if !isnothing(window_of_interest)
                    pre::Int  = window_of_interest[1] 
                    post::Int = window_of_interest[2] 
                    signal = signal_matrix[ch_ind, si-pre:si+post]
                    ind_start = si-pre
                else
                    signal = signal_matrix[ch_ind, cum_len:cum_len+fl-1]
                    ind_start = cum_len
                end
                
                pk_i, pk_prop = findpeaks1d(signal,
                    height = pk_bkgnds[fn][ch_ind],
                    distance = intra_chan_window)

                pk_i .+= ind_start
                largest_peaks[pk_i] = max.(largest_peaks[pk_i], pk_prop["peak_heights"])
            end
        end
        cum_len += fl
    end
    
    pk_inds, pk_hts = findpeaks1d(largest_peaks, distance = inter_chan_window, height = 0.0)
    pk_hts =  pk_hts["peak_heights"] 

    #? AF edited
    #Added line to recreate pk_inds based on waveforms within the stimulus window 
    #such that length(pk_inds)=size(extracted_peaks,2)
    pk_inds_new = Vector()  #? AF added to create new vector
    extracted_peaks = Vector()
    for pk_ind in pk_inds
        pre = pk_ind - extract_window[1]
        post = pk_ind + extract_window[2]
        if (pre <= 0 || post > sig_len)
            continue
        end
        pk_sample = vcat([signal_matrix[chi,pre:post] for chi = 1:num_channels]...)
        push!(extracted_peaks, pk_sample)
        push!(pk_inds_new,pk_ind)   #?AF: Added line to push only the pk_inds associated with extracted peaks

    end

    concat_sigs = hcat(extracted_peaks...)
    return pk_inds_new, concat_sigs, pk_hts
end

function combine_templates(retemplate, assignsOriginal, templatesOriginal, assignsNew, templatesNew, templatecounter, plot_output = true)
    # Order all assignments according to count size
    # assignsOrder is a vector of pairs (bool for new or old (0 = old, 1 = new) and int for assignment number)
    # Initialise assignsOrder with original assigments larger than the retemplated assignment, as they must have the greatest overall counts
    assignsOrder = [(false => x) for x in 1:(retemplate - 1)]
    n = 1
    # Iterate over remaining original assignments (excluding retemplated)
    for originalAssign in (retemplate + 1):length(templatesOriginal[1,:])
        # While next new assignment has greater count than next original assignment, push the new assignment
        # First check that there are new assignments remaining (i.e. not trying to access counts(assignsNew) at too large an index)
        while n <= length(counts(assignsNew)) && counts(assignsNew)[n] > counts(assignsOriginal)[originalAssign]
            push!(assignsOrder, (true => n))
            n += 1
        end
        # Otherwise push next original assignment
        push!(assignsOrder, (false => originalAssign))
    end
    # Push any remaining new assignments
    for newAssign in n:length(templatesNew[1,:])
        push!(assignsOrder, (true => newAssign))
    end

    # Combine templates
    templatesVec = []                                                   # construct empty vector
    for assign in assignsOrder
        if assign.first                                                 # if assignment from assignsNew
            push!(templatesVec, templatesNew[:, assign.second])         # push from templatesNew
        else push!(templatesVec, templatesOriginal[:, assign.second])   # else push from templatesOriginal
        end
    end
    templates = reduce(hcat, templatesVec)

    # Combine assignments
    n = 1                                                               # clear n
    assigns = Int64[]                                                   # construct empty vector
    templatecounternew = zeros(length(unique(assignsOriginal)) + length(unique(assignsNew)) -1) #make a new template counter and fill it with zeroes

    for assign in assignsOriginal
        # If assign was retemplated, push from assignsNew
        if assign == retemplate
            # Check that there is only one index for the relevant assign in assignsOrder
            assignmentIndex = findall(x -> x == (true => assignsNew[n]), assignsOrder)
            if length(assignmentIndex) == 1
                push!(assigns, assignmentIndex[1])
                templatecounternew[assignmentIndex[1]] = templatecounter[assign] .+ 1 #iterate the template counter to signal that these assignments have been retemplated
            # If there is not exactly 1 index, throw an error
            else
                repeatedAssign = assignsNew[n]
                error("The assignsOrder vector contains duplicates of (1 => $repeatedAssign)")
            end
            n += 1
        # Else push from assignsOriginal (same steps as above)
        else
            assignmentIndex = findall(x -> x == (false => assign), assignsOrder)
            if length(assignmentIndex) == 1
                push!(assigns, assignmentIndex[1])
                templatecounternew[assignmentIndex[1]] = templatecounter[assign] # don't interate the template counter
            else
                error("The assignsOrder vector contains duplicates of (0 => $assign)")
            end
        end
    end
    println("Combined counts: ", counts(assigns))

    # Generate plots for each new assignment to check the output
    
    if plot_output
        for waveNo in 1:length(templates[1,:])
            temporaryWave = fndAssign(waveforms, waveNo, assigns)
            figure()
            PyPlot.plot(temporaryWave, color=COLOURS[waveNo])
            PyPlot.plot(templates[:, waveNo], color="black")
        end
    end
    return assigns, templates, templatecounternew
end

function isiTests(spikes, onems; reject=true, trials_threshold=0.5, isi_threshold=1, percentage_threshold=1)
    # Spike number test
    spikes_count = sum(length(x) for x in spikes)
    trials_count = length(spikes)
    spike_test_fail = spikes_count < trials_threshold * trials_count
    println("$spikes_count spikes in $trials_count trials")
    
    # ISI test
    isi_fails = []
    for trial in spikes
        for (i, spike) in enumerate(trial)
            if i < length(trial) && trial[i+1] - spike <= isi_threshold*onems
                push!(isi_fails, spike)
            end
        end
    end
    isi_fails_count = length(isi_fails)
    percentage_isi_failed = round(isi_fails_count/spikes_count*100, digits = 2)
    isi_test_fail = percentage_isi_failed >= percentage_threshold
    println("$isi_fails_count ISIs <= $isi_threshold ms ($percentage_isi_failed%)")

    # Return test results
    if spike_test_fail && reject
        return "REJECTED"
    elseif isi_test_fail
        return "FAILED"
    else return "PASSED"
    end
end

"""
Takes a row from the sOutput spreadsheet (in DataFrame format) and a string (either "Ref" or "Text") and returns non-global variables for an input file.
"""
function loadSOutputFile(df_row::DataFrameRow{DataFrame, DataFrames.Index}, refOrTest::String)
    filename = df_row[refOrTest*"Filename"]                                                     # name of file to be loaded, used in some error messages
    file = addExtToFilename ?                                                                   # check if file extension should be appended
            load("spikesortingOutput/$filename.jld2") : load("spikesortingOutput/$filename")    # loads the relevant file
    pksLocs = file["pksLocs"]
    pksAss = file["pksAss"]
    ass_counts = file["ass_counts"]
    pre_trial_len = file["pre_trial_len"]                                                       # time before stimulus
    post_trial_len = file["post_trial_len"]                                                     # time after stimulus
    pks_plot = file["pks_plot"]                                 
    onems = file["onems"]                                                                       # samples per millisecond
    stim_code = file["stim_type"]                                                               # code for stim type
    bad_assigns = file["badAssigns"]                                                           # list of assignments to be ignored
    return filename, file, pksLocs, pksAss, ass_counts, pre_trial_len, post_trial_len, pks_plot, onems, stim_code, bad_assigns
end

"""
    function reorder_kcr(kcr)

    Reorders the assignment in the Kmeans clustering result (kcr) so that the largest group of values are assigned as 1, second largest 2, and so on

    The sp provides the assignment order for modifying template outout in cluster spikes (via apply_kmeans)
"""

function reorder_kcr(kcr)
    new_assigns = zeros(Int,length(kcr.assignments)) #Making a copy of the cluster assignments
    new_centers = zeros(size(kcr.centers))  #Making a copy of the centers
    assign_order = sortperm(counts(kcr.assignments),rev=true)  #sortperm is the order to sort the values in, we want rev so that its largest first
    for (new,original) in enumerate(assign_order)
        new_assigns[kcr.assignments .== original] .= new
        new_centers[:,new] = kcr.centers[:,original]
    end
    return KmeansResult(new_centers, new_assigns, kcr.costs, kcr.counts, kcr.wcounts, kcr.totalcost, kcr.iterations, kcr.converged), assign_order 
end