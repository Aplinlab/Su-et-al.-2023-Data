Lbv = v"4.2.2+TS"
# Batch version of loading.
#! THIS VERSION ONLY WORKS IF YOU ALREADY KNOW WHAT PARAMETERS TO FEED INTO LOADING USING THE LOOKUP SPREADSHEET

# * Changelog
#   Completely different for loop
#   Formatting/order changes
#   More error messages
#   Commented out some obsolete/plotting-only variables
#   Define trials_per_file for laser according to files_num
#   See loading for additional changelog

# ! Errors

#   TODO List
#   Move regex out of loop
#   Clean up for upload

# ? Suggestions & Queries



# ! before loading new files, ensure you kill the julia session to remove all global variables
# Packages
using DataFrames
using Query
using XLSX
using PyPlot; pygui(true)
using Statistics # std
using StatsBase # mode
using JLD2

# Working directory and dependencies
cd(dirname(@__FILE__))
include("src/IntanReader.jl")
include("src/preprocessing_module_TS1.1.jl")        # Preprocessing Tools
include("src/response_quantification_TS1.5.jl")     # Response quantification Tools
include("src/RHS2000_MEA_config_v1.jl")



df_key = DataFrame(XLSX.readtable("savedFilesInfo2023v3.xlsx", "lOutput and sTemplate"))
df_inverter = DataFrame(XLSX.readtable("savedFilesInfo2023v3.xlsx", "Descriptions"))

# filename_spreadsheet = "metadataJackFPA.xlsm"
filename_spreadsheet = "metadataTomTFSv4.xlsm"
#datafile_location    = "R:/Raw Data/2021 iDC Data/"     # Tom and Felix uni
# datafile_location    = "R:/Potas/Raw Data/2021 iDC Data/"   # Jason uni
datafile_location    = "Raw Data/"   # Tom home

logic_high = true               # whether the triggers are peaks or valleys
plot_all = true                 # false will produce a single waveform plot (for use as a figure)

# Stim Artefact
art_pre = 60            # use to control pre-stim blanking period (30 = 1 ms)
art_post = 60           # use to control post-stim blanking period
art_height = 500        # used to threshold artefact blanking (when using artefact blank only)
trig_blank = true       # To control whether the artefact itself is used to blank, or the digital trigger timing (true = trigger)

# Filtering
bpass = 300
bstop = 5000

# Peak Detection
intra_chan_window = 5
inter_chan_window = 6
peak_assignment_err = 5 # 0,3,5
discard_first = 100                 # NOTE Remove first X samples from signal due to recording artefact
elec_stim_thresh = 500 # 180        #Threshold to find peaks
mechan_stim_thresh = 0.75

# Fidelity Computation (to find Best Responses, i.e. reproducibility):
# kernel_size = 0.3 / 1000            # fs
# trial_ratio = 0.4
# bin_per_ms = 1                      # NEW #1 means `responses` will have 1200 length
# fidelity_matrix_reorder_channels = true

for loadinginput in eachrow(df_key)

    output_filename = loadinginput["Filename"]
    template_filename = loadinginput["Template"]

    animalID = loadinginput["Animal"]
    pos_sheet = loadinginput["Position"]
    cell_sheet = loadinginput["Cell"]

    df_inverter_t = df_inverter |>
        @filter(_.Template == template_filename) |>
        @filter(_.Animal == animalID) |>
        @filter(_.Position == pos_sheet) |>
        @filter(_.Cell == cell_sheet) |>
    DataFrame
    # Throw error if more than one matching template is found
    if size(df_inverter_t)[1] > 1
        error("More than one template found in Descriptions for $animalID $output_filename")
    elseif size(df_inverter_t)[1] == 0
        error("No matching template found in Descriptions for $animalID $output_filename")
    end
    invert_bool = df_inverter_t[1, "trace_inverter"]
    if invert_bool == "y"   #set up inverter logic
        trace_inverter = ["all"]
    elseif invert_bool == "n" 
        trace_inverter = []
    else error("Invalid trace_inverter value for $animalID $output_filename")
    end

    # * Parameters
    # Regex
    regexParam = r"(?P<position>\d+) (?P<amp>\d\.?\d*) (?P<block_amp>-?\d+) (?P<stim_type>\w{2}) (?P<files_num>\[(?:\d+(?:, )*)+]) (?P<channels>\[(?:\d+(?:, )*)+])"
    parameters_output = match(regexParam, output_filename)
    parameters_template = match(regexParam, template_filename)

    # Error checking
    pos_template = String(parameters_template["position"])
    pos_output = String(parameters_output["position"])
    if pos_sheet != pos_template || pos_sheet != pos_output
        error("Positions do not match for $animalID $output_filename")
    end
    stim_type_sheet = df_inverter_t[1, "Stimtype"]
    stim_type_template = String(parameters_template["stim_type"])
    if stim_type_sheet != stim_type_template
        error("Template stim types do not match for $animalID $output_filename")
    end
    channels_sheet = df_inverter_t[1, "Channels"]
    channels_template = String(chop(parameters_template["channels"]; head=1, tail=1))
    channels_output = String(chop(parameters_output["channels"]; head=1, tail=1))
    if channels_sheet != channels_template || channels_sheet != channels_output
        error("Channels lists do not match for $animalID $output_filename")
    end

    # Define Parameters
    pos = pos_sheet
    amp = String(parameters_output["amp"])
    block_amp = String(parameters_output["block_amp"])
    stim_type = String(parameters_output["stim_type"])
    files_num = parse.(Int, split(chop(parameters_output["files_num"]; head=1, tail=1), ','))

    # IPI = "*" #! note filter for IPI is off below
    # SpC_side = "L"
    # limb = "LH"
    # channels = [35, 106, 108]
    # WOI = [10,100]

    # Additional Parameters
    if stim_type == "L1"
        pretrial_time = 1.0
        trial_t = 5.0
        pre_ratio_window = [-500, 0]
        post_ratio_window = [1000, 1500]
        binwidth = 50
        bkg_thresh = 3.0
        if length(files_num) > 2
            error("Too many files_num values for a laser set ($animalID $output_filename)")
        end
        global trials_per_file = 3 - length(files_num)
    elseif stim_type == "M3"
        pretrial_time = 0.5
        trial_t = 2.0
        pre_ratio_window = [-500, 0]
        post_ratio_window = [50, 550]
        binwidth = 10
        bkg_thresh = 3.0
        global trials_per_file = 10
    elseif stim_type == "M2"
        pretrial_time = 0.05
        trial_t = 0.4
        pre_ratio_window = [-1, 0]
        post_ratio_window = [0, 1]
        binwidth = 1
        bkg_thresh = 3.5
        global trials_per_file = 10
    else
        pretrial_time = 0.02
        trial_t = 0.1
        pre_ratio_window = [-1, 0]
        post_ratio_window = [0, 1]
        binwidth = 1
        bkg_thresh = 3.5
        global trials_per_file = 10
    end

    # Response Quantification
    # resp_ratio_thresh = 2               # Ratio of peaks after stimulus to before stimulus #What does resp mean?
    # resp_diff_thresh = 7 * length(files_num)    # Difference of peaks after stimulus compared to before stimulus

    ################################################################################

    # Filtering data
    # ' ## Getting filenames
    # Spreadsheet should have animalID as sheet names
    df = DataFrame(XLSX.readtable(datafile_location*filename_spreadsheet, animalID, stop_in_empty_row=false))
    df_t = df |>
        @filter(_.animal == animalID) |>
        # @filter(_.SpC_side == SpC_side) |>
        # @filter(_.timing == IPI) |>
        # @filter(_.limb == limb) |>
        @filter(_.electrode_pos == pos) |>
        @filter(_.amplitude == amp) |>
        @filter(_.block_amp == block_amp) |>
        @filter(_.stim_type == stim_type) |>
        @map({_.rhs_name, _.electrodes}) |>
    DataFrame
    # NOTE: make sure your spreadsheet is correctly formatted (e.g. no spaces after )

    ################################################################################

    # + echo=false; results="hidden"

    # include("D:\\Documents\\OneDrive - UNSW\\1-research\\active projects\\Lab2021\\Spike sorting package\\src\\IntanReader.jl")
    # include("D:\\Documents\\OneDrive - UNSW\\1-research\\active projects\\Lab2021\\Spike sorting package\\src\\preprocessing_module_v1.jl")  # Preprocessing Tools
    # include("D:\\Documents\\OneDrive - UNSW\\1-research\\active projects\\Lab2021\\Spike sorting package\\src\\response_quantification_v1.jl") # Response quantification Tools

    # include rhsMEA_config                 # Probe configuration and Channel ordering

    probes = df_t[!,:electrodes]
    rhsMEA_config(probes[1])

    filenames = convert(Array, df_t[:,1][files_num]) .* ".rhs"
    @info filenames  # tells what is in filenames
    probes = union(df_t[!,:electrodes]) # checks to ensure all probes are same for the filtered data set
    if length(probes) == 1
        probe_config = first(probes)
        probes = first(probes)
    else
        @error "Selected files have different probe arrangement"
        probes = first(probes)
    end

    active_chs = sort(filter(!iszero, rhsMEA_config(probes)[:]))
    # active_chs = sort(filter(!iszero, rhsMEA_config(probes[1])[:]))

    max_channels = maximum(active_chs) # How many channels have recordings

    # ' ## Loading Data
    # Variable Setup
    sig_neg_vec = []
    stim_command_vec = []
    stim_timing_vec = []
    laser_vec = []
    file_lengths = []
    bad_channels = []

    # fn = filenames[1]
    for fn in filenames
        rhs_name = datafile_location * animalID * "/" * fn
        read_data(rhs_name)

        # Labelling Loaded Data
        file_length = length(amplifier_data[1, discard_first:end])
        sig_neg = zeros(max_channels, file_length)
        sig_neg[active_chs,:] = -amplifier_data[active_chs, discard_first:end]
        # laser = board_adc_data[1, discard_first:end]

        if stim_type == "M3" # if using piezo stimulation, use piezo stim timing function and don't remove artefacts
            trigthresh = 0.1
            pzChannel = 1
            stim_timing = []
            stim_command = []
            pz_trig = board_adc_data[pzChannel, :]

            filteredTrigsig = filt(digitalfilter(Lowpass(5000; fs = 30000), Butterworth(4)), pz_trig)
            for i in 1:100
                avtrig = mean(filteredTrigsig[10:20000])
                filteredTrigsig = filteredTrigsig .- avtrig
            end
        
            for sampleNo in 10001:length(pz_trig)
                if filteredTrigsig[sampleNo] > trigthresh && maximum(filteredTrigsig[(sampleNo - 10000):(sampleNo - 1)]) <= trigthresh
                    push!(stim_timing, sampleNo-3000)
                    push!(stim_command, sampleNo-3000)
                end
            end
        else
            stim_command = logic_high ? board_dig_in_data[1, discard_first:end] : 1 .- board_dig_in_data[1, discard_first:end] # Not used
            stim_timing = logic_high ? board_dig_in_data[2, discard_first:end] : 1 .- board_dig_in_data[2, discard_first:end]
        end

        # Can close file here
        push!(sig_neg_vec, sig_neg)
        push!(stim_command_vec, stim_command)
        push!(stim_timing_vec, stim_timing)
        # push!(laser_vec, laser)
        push!(file_lengths, file_length)
    end

    global fs = frequency_parameters.amplifier_sample_rate
    global post_trial_t = trial_t - pretrial_time # post-stim trial time (s)
    global pre_trial_len = convert(Int, fs * pretrial_time)
    global trial_len = convert(Int, fs * trial_t)
    global post_trial_len = convert(Int, round(fs * post_trial_t))
    global t_trial = (1000 * (1:trial_len) ./ fs) .- 300 # times in ms for plotting

    global total_file_length = sum(file_lengths)
    global num_files = length(file_lengths)
    # global half_height = minimum([num_files*trials_per_file, 10])
    println("Loading finished")

    # ' ## Preprocessing Data

    # 1) Finding Stimulus Timing
    stim_inds_vec = []
        for fno = 1:num_files
            if stim_type == "M3" # if using piezo stimulation, use piezo stim timing function and don't remove artefacts
                stim_inds = [stim_timing_vec[fno] for ch = 1:max_channels]
                new_bad_channels = []

            # elseif stim_type[1] == 'E' # Electrical stimulus #! Deprecated code for triggering off the electrical artefact
            #         # Find stim inds accepts DxN matrix, so sig1_neg needs to be transposed
            #         # stim_inds, new_bad_channels = find_stim_inds(sig_neg_vec[fno]', elec_stim_thresh) # this is for triggering off the singal
            #         sis, new_bad_channels = find_stim_inds(stim_timing_vec[fno], mechan_stim_thresh) # this is for triggering off the pulse trigger
            #         stim_inds = [sis[1] for ch = 1:max_channels]
        
            else  # otherwise, use the stim timing from the sync trigger and remove artefacts
                sis, new_bad_channels = find_stim_inds(stim_timing_vec[fno], mechan_stim_thresh)
                stim_inds = [sis[1] for ch = 1:max_channels]
            end
            
            global active_chs = filter(x -> x ∉ new_bad_channels, active_chs)
            push!(stim_inds_vec, stim_inds)
        end
    # stim_inds_vec = []
    # for fno = 1:num_files
    #     if stim_type[1] == 'E' # Electrical stimulus
    #             # Find stim inds accepts DxN matrix, so sig1_neg needs to be transposed
    #             # stim_inds, new_bad_channels = find_stim_inds(sig_neg_vec[fno]', elec_stim_thresh) # this is for triggering off the singal
    #             sis, new_bad_channels = find_stim_inds(stim_timing_vec[fno], mechan_stim_thresh) # this is for triggering off the pulse trigger
    #             stim_inds = [sis[1] for ch = 1:max_channels]
    
    #         else
    #             sis, new_bad_channels = find_stim_inds(stim_timing_vec[fno], mechan_stim_thresh)
    #             stim_inds = [sis[1] for ch = 1:max_channels]
    #         end
    #         global active_chs = filter(x -> x ∉ new_bad_channels, active_chs)
    #         push!(stim_inds_vec, stim_inds)
    #     end

    # 2) Remove Artefacts where applicable (i.e. stim_type == "E" or if you have a sync signal)
    sig_rm_art_vec = []
    for fno = 1:num_files # file number
        sig_rm_arts = deepcopy(sig_neg_vec[fno])
        sig_rm_trigs = deepcopy(stim_inds_vec[fno])
        if trig_blank == false #remove artifacts using threshold   
                if stim_type[1] == 'E'  # artefact will be removed from E1 and M1
                # Using remove_artefacts() from Preprocessing Tools
                sig_rm_arts[:,pre_trial_len:end], new_bad_channels = remove_artefacts(sig_rm_arts[:,pre_trial_len:end])
                append!(bad_channels, new_bad_channels)
            end

        elseif trig_blank == true #remove artifacts using the trigger signal
            
            #define some variables for the interpolation
            art_len = art_pre + art_post-1
            lin_range = LinRange(1,2, art_len)
            num_chan = size(sig_rm_arts,1)
            
            for ch = 1:num_chan
                art_ind = sig_rm_trigs[ch] #find the trigger timings
                for i in art_ind #where each i is an artefact index
                    start_pos = i - art_pre
                    interpol = interpolate([sig_rm_arts[ch,start_pos],sig_rm_arts[ch,start_pos+art_len]], BSpline(Quadratic(Periodic(OnGrid())))) #interpolate out the area around the sync trigger
                        for r in 1:art_len #where each r is the sample index in the artefact duration
                            sig_rm_arts[ch,start_pos+r-1] = interpol(lin_range[r])
                        end
                end
            end
            

        end    
        # If not electrical stim, do nothing

        # Inverter code to selectively invert channels
        if trace_inverter == []
        #If no inverter values, skip code
        elseif trace_inverter == ["all"]
        #If inverter is "on", will invert all channels
            sig_rm_arts = -sig_rm_arts
        else                        
            for n = 1:length(trace_inverter)
            #The main for loop to go through the selected channels to invert    
                if trace_inverter[n] > 128 || trace_inverter[n] < 1
                #If statement to ignore input errors, statement selective to a 128 electrode array
                    l = trace_inverter[n]
                    println("$l is not a valid channel")
                else                    
                    sig_rm_arts[trace_inverter[n],:] = -sig_rm_arts[trace_inverter[n],:]
                    #Converting channel to negative
                end    
            end
        end
        push!(sig_rm_art_vec, sig_rm_arts)
    end
    # 3) Filtering
    HFdata_neg_vec = []
    for fno = 1:num_files
        # Using filter_signal() from Preprocessing Tools
        HFdata_neg = filter_signal(sig_rm_art_vec[fno])
        push!(HFdata_neg_vec, HFdata_neg)
    end
    global active_chs = filter(x -> x ∉ bad_channels, active_chs) # filter for active channels that aren't in the bad set

    #=     # 3) Finding Stimulus Timing
    stim_inds_vec = []
    for fno = 1:num_files
        if stim_type[1] == 'E' # Electrical stimulus
                # Find stim inds accepts DxN matrix, so sig1_neg needs to be transposed
                # stim_inds, new_bad_channels = find_stim_inds(sig_neg_vec[fno]', elec_stim_thresh) # this is for triggering off the singal
                sis, new_bad_channels = find_stim_inds(stim_timing_vec[fno], mechan_stim_thresh) # this is for triggering off the pulse trigger
                stim_inds = [sis[1] for ch = 1:max_channels]
    
            else
                sis, new_bad_channels = find_stim_inds(stim_timing_vec[fno], mechan_stim_thresh)
                stim_inds = [sis[1] for ch = 1:max_channels]
            end
            global active_chs = filter(x -> x ∉ new_bad_channels, active_chs)
            push!(stim_inds_vec, stim_inds)
        end
    =#
    # 4) Defining background activity for each channel
    bkgnd_means = []
    bkgnd_stds = []
    peak_backgrounds = []
    for fno = 1:num_files
        bkgnd_mean = zeros(max_channels)
        bkgnd_std = zeros(max_channels)
        for ch in active_chs
            stim_ind = stim_inds_vec[fno][ch]
            be::Int = first(stim_ind) - 100 # background end
            if be < 1 * fs
                @warn "Less than 1s of background available"
            end
            bs::Int = max(1, be - fs) # background start
            bm = mean(HFdata_neg_vec[fno][ch, bs:be])
            bstd = std(HFdata_neg_vec[fno][ch, bs:be])
            bkgnd_mean[ch] = bm
            bkgnd_std[ch] = bstd
        end
        push!(bkgnd_means, bkgnd_mean)
        push!(bkgnd_stds, bkgnd_std)
        push!(peak_backgrounds, bkgnd_mean + (bkg_thresh * bkgnd_std)) # Select peaks with 3 s.d. above mean
    end

    # ' ## Split into trials
    ####################### Split into Trials ##############################
    # ' This allows graphing trials separately
    # Ch_trials_data = [Vector() for ch in 1:max_channels]
    # Ch_trials_spikes = [Vector() for ch in 1:max_channels]
    # Ch_trials_spikes_h = [Vector() for ch in 1:max_channels]
    # for fno = 1:num_files
    #     for (i, ch) in enumerate(active_chs)
    #         pk_thresh = bkgnd_stds[fno][ch] * bkg_thresh # 3 times standard deviation
    #         # Using split_trials() from Preprocessing Tools
    #         trial_data = split_trials(HFdata_neg_vec[fno][ch,:], stim_inds_vec[fno][ch])
    #         trial_spikes = Vector()
    #         trial_spikes_h = Vector()
    #         for trial in trial_data
    #             spikes, spike_prop = scisig.find_peaks(trial, height=pk_thresh, distance=intra_chan_window)
    #             push!(trial_spikes, spikes)
    #             push!(trial_spikes_h, spike_prop["peak_heights"])
    #         end
    #         # Ch_trial_data[ch] = trial_data
    #         append!(Ch_trials_data[ch], trial_data)
    #         # Ch_trial_spikes[ch] = trial_spikes
    #         append!(Ch_trials_spikes[ch], trial_spikes)
    #         # Ch_trial_spikes_h[ch]
    #         append!(Ch_trials_spikes_h[ch], trial_spikes_h)
    #     end
    # end

    ####################### Response Quantification ##############################
    # response_ratio_Mat = fill(NaN, 128)
    # response_diff_Mat = fill(NaN, 128)
    # response_spikeFreq_Mat = [Vector() for ch in 1:max_channels]
    # total_num_trials = trials_per_file * num_files
    # # num_channels = length(active_chs)

    # bin_edges = collect(-pretrial_time * 1000:binwidth:(-pretrial_time + trial_t) * 1000)
    # trial_len = round(Int, fs * trial_t)
    # t_trial = (1000 * (1:trial_len) ./ fs) .- 300
    # fidelity_base = zeros(max_channels, Int(fs * trial_t))

    # for ch in active_chs
    #     Ch_trial_spikes = Ch_trials_spikes[ch]
    #     for peaks in Ch_trial_spikes
    #         fidelity_base[ch, peaks] .+= 1
    #     end
    #     spikemat_concat = vcat(Ch_trial_spikes...) ./ (fs / 1000)
    #     spike_freq, spike_binEdges = np.histogram(spikemat_concat .- 300, bins=bin_edges)
    #     resp_ratio, resp_diff = response_quantify(spike_freq, spike_binEdges, pre_ratio_window, post_ratio_window)
    #     response_ratio_Mat[ch] = resp_ratio
    #     response_diff_Mat[ch] = resp_diff
    #     response_spikeFreq_Mat[ch] = spike_freq
    # end

    # # Fidelity Computation
    # kernel = ones(round(Int, kernel_size * fs))
    # convul = conv(kernel, fidelity_base')'
    # convul_norm = convul ./ total_num_trials
    # fidelity = zeros(max_channels, Int(trial_t * 1000 * bin_per_ms) + 1)
    # for (i, ch) in enumerate(active_chs)
    #     c_sig, _ = scisig.find_peaks(convul_norm[i,:], distance=intra_chan_window)
    #     c_sig_real = c_sig[ (c_sig .> fs / (1000 * bin_per_ms)) .& (convul_norm[i,c_sig] .> trial_ratio)]
    #     fidelity[ch,round.(Int, c_sig_real ./ fs .* 1000 .* bin_per_ms)] .= convul_norm[i,c_sig_real]
    # end

    # ############################ Showing Response Quantification ######################
    # # Fidelity Matrix with train data
    # pyplot_fidelity_matrix(fidelity, reorder = false) # reorder can't work without RHSMEA 2D probe configuration

    # laser_Mat_concat_plot = setup_laser_data(laser_vec, stim_inds_vec)
    # Ch_trials_data_plot = Ch_trials_data
    # Ch_trials_spikes_plot = Ch_trials_spikes
    # Ch_trials_spikes_h_plot = Ch_trials_spikes_h
    # response_ratio_Mat_plot = response_ratio_Mat
    # response_diff_Mat_plot = response_diff_Mat
    # pyplot_ratio_diff_grid(response_ratio_Mat[1:max_channels], response_diff_Mat[1:max_channels])

    ################################################################################

    println("Loading Complete")


    channels = parse.(Int, split(chop(parameters_output["channels"]; head=1, tail=1), ','))
    signal_matrix = hcat(HFdata_neg_vec...)[channels,:]
    pk_bkgnds = [ file[channels] for file in peak_backgrounds]
    file_lengths
    cum_file_length = cumsum(file_lengths)

    stim_inds = [stimulus_indices[channels[1]] for stimulus_indices in stim_inds_vec]

    cum_stim_inds = deepcopy(stim_inds_vec[1][channels[1]])
    for fv in 2:num_files
        append!(cum_stim_inds, copy(stim_inds_vec[fv][channels[1]]) .+ cum_file_length[fv - 1])
    end

    mkpath("loadingOutput/$animalID")
    # @save "loadingOutput/$animalID/$output_filename.jld2" signal_matrix pk_bkgnds file_lengths stim_inds cum_stim_inds animalID SpC_side limb IPI pos amp block_amp stim_type filenames files_num channels
    @save "loadingOutput/$animalID/$output_filename.jld2" signal_matrix pk_bkgnds file_lengths stim_inds cum_stim_inds animalID pos amp block_amp stim_type filenames files_num channels

    println("Loading Saved")
end