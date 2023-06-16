"""
    preprocessing_module_v1.jl
    see loading for changelog

    preprocess and feature selection

    * filter_sgnal(signal) #NxD
    * function remove_artefacts(signal; debug=true) #NxD
    * find_stim_inds(signals, stim_thresh) #DxN
    * split_trials(signal, stim_inds) #D
    * assign_to_trials(stim_inds, peak_inds, iterable) #N

    Leave peak detection top level to understand indices - trials or array
    * Peak
"""

using DSP
using Interpolations
using PyCall
scisig = pyimport("scipy.signal")

"""
    function filter_signal(signal)
    global bpass, bstop, fs

    Filters the signal with a bandpass digital filter
    `signal` NxD matrix. N signals of D length
    Returns a NxD matrix.
"""
function filter_signal(signal)
    global bpass, bstop, fs
    responsetype = Bandpass(bpass, bstop, fs=fs)
    designmethod = Butterworth(3)
    filtered_signal = filtfilt(digitalfilter(responsetype, designmethod), signal')
    return filtered_signal'
end

"""
    function remove_artefacts(signal; debug=true)

    Uses global bad_channels, art_pre, art_post, art_height, trial_t, fs
    Removes artefacts in the signal which are above `art_height`
    `signal` is of size(N,D). N channels of D samples
    Returns `signal` (changes made inplace)
"""
function remove_artefacts(signal; debug=true)
    global bad_channels, art_pre, art_post, art_height, trial_t, fs
    #should I check bad channels?
    new_bad_channels = Vector()
    min_dist = trial_t * fs /2
    num_chan = size(signal,1)
    sig_len = size(signal,2)
    art_len = art_pre + art_post-1
    for ch = 1:num_chan
        art_ind, _ = scisig.find_peaks(signal[ch,:], height = art_height, distance = min_dist)
        if length(art_ind) == 0
            push!(new_bad_channels, ch)
            if debug
                @warn "Remove Artefact: No artefact with height $art_height on channel $ch"
            end
        elseif art_ind[1] - art_pre < 1
            push!(new_bad_channels, ch)
            if debug
                @warn "Remove Artefact: First artefact found in first $art_pre samples of signal"
            end
        elseif art_ind[end] + art_post > sig_len
            push!(new_bad_channels, ch)
            if debug
                @warn "Remove Artefact: Last artefact found in last $art_post samples of signal"
            end
        else
            lin_range = LinRange(1,2, art_len)
            for i in art_ind #where each i is an artefact index
                start_pos = i - art_pre
                interpol = interpolate([signal[ch,start_pos],signal[ch,start_pos+art_len]], BSpline(Quadratic(Periodic(OnGrid()))))
                for r in 1:art_len #where each r is the sample index in the artefact duration
                    signal[ch,start_pos+r-1] = interpol(lin_range[r])
                end
            end
        end
    end
    # append!(bad_channels, new_bad_channels)
    return signal, new_bad_channels
end

# function channels_find_peaks(signal, num_channels, thresh = zeros(num_channels))
#     return thresh
# end

"""
    function find_stim_inds(signals, stim_thresh)
    Uses global trials_per_file to determine how many stim_inds there are per channel

    `signals` is DxN matrix or D vector
    Returns Vector of length D, of vectors with `trials_per_file` length
"""
function find_stim_inds(signals, stim_thresh)
    global pre_trial_len, post_trial_len, trials_per_file
    local bad_channels = Vector()
    stim_inds = Vector()
    ch = 1
    for signal in eachcol(signals)
        pk, pkprop = scisig.find_peaks(signal[pre_trial_len: (end-post_trial_len)], height = stim_thresh, distance = (pre_trial_len+post_trial_len)/2, plateau_size = 10)
        stim_ind = pkprop["left_edges"]
        stim_ind .+= pre_trial_len
        stim_ind_real = stim_ind[sortperm(pkprop["peak_heights"],rev=true)][1:min(end,trials_per_file)]
        push!(stim_inds, sort(stim_ind_real))
        if length(stim_ind_real) != trials_per_file
            push!(bad_channels, ch)
        end
        ch += 1
    end
    return stim_inds, bad_channels
end

"""
    function split_trials(signal, stim_inds)

    Uses global variables: fs, pre_trial_len, post_trial_len
    Splits the signals into trials based off indices in `stim_inds`
"""
function split_trials(signal, stim_inds)
    global fs, pre_trial_len, post_trial_len
    trials = Vector()
    for stim_ind in stim_inds
        sample = signal[(stim_ind-pre_trial_len):(stim_ind+post_trial_len-1)]
        push!(trials, sample)
    end
    return trials
end

"""
    function assign_to_trials(stim_inds, peak_inds, iterable)

    Uses global pre_trial_len, post_trial_len
    The values of peak_inds should be in the same scale as with stim_inds in a single vector
    Each stim_ind becomes a "trial" and the peak_inds are separated into these trials
    The number of elements in `iterable` should be the same as in `peak_inds`
    USAGE EXAMPLES:
    - assign_to_trials(stim_inds, peak_inds, peak_heights)
    - assign_to_trials(stim_inds, peak_inds, peak_assigns)
"""
function assign_to_trials(stim_inds, peak_inds, iterable)
    global pre_trial_len, post_trial_len
    peak_ind_trials = Vector()
    peak_ind_properties = Vector()
    for stim_ind in stim_inds
        peak_ind_trial = Vector()
        peak_ind_prop = Vector()
        trial_start = (stim_ind-pre_trial_len)
        trial_end = stim_ind+post_trial_len-1
        for (i,peak_ind) in enumerate(peak_inds)
            if peak_ind > trial_end
                break #Move to next trial/stim_ind
            end
            if peak_ind in trial_start:trial_end
                # Reset value
                push!(peak_ind_trial, peak_ind - trial_start +1)
                push!(peak_ind_prop, iterable[i])
            end
        end
        push!(peak_ind_trials, peak_ind_trial)
        push!(peak_ind_properties, peak_ind_prop)
    end
    return peak_ind_trials, peak_ind_properties
end
