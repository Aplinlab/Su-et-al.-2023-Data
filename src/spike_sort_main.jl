"""
    Spike sorting functions for template generation and template matching
"""


COLOURS = ["green","blue","red", "y", "orange", "k", "pink", "gray",
    "olive", "orange", "#8c564b", "#e377c2"]

using Clustering # for CLustering
using Distances # for Clustering
using PyCall # for Spike detection
using PyPlot # for Display
using Printf # for Display
using MultivariateStats # for PCA
using Random # for Clustering
using FindPeaks1D # for native julia peakfinding
using StatsBase #mode
# scisig = pyimport("scipy.signal")
# MultiCursor = matplotlib.widgets.MultiCursor

"""
    function detect_spikes(signal_matrix::AbstractMatrix,
        pk_bkgnds::AbstractVector, file_lengths::AbstractVector, stim_inds, window_of_interest;
        fs=30000,
        intra_chan_window=6,
        inter_chan_window=6,
        extract_window = [10,10])

    Detects spikes in the signal_matrix and returns
    * The index locations and heights of the spikes
    * The matrix of extracted spike waveforms

    Arguments
    `signal_matrix` NxD matrix of the signal where the rows are the channels (N)
    `pk_bkgnds` A vector of thresholds for detecting a peak, each element should be a vector
                of length N
    `extract_window` 2 element array indicating window around spike for extraction. Default is [10,10]
                (10 before and 10 after the spike, total 21 samples) (0.7ms)

    `file_lengths` A vector of lengths which help break down the signal_matrix and apply the pk_bkgnd
                at the correct interval. These intervals are assumed to have come from different "files"
    `stim_inds` The stimulus indices from each "file".
                # If not provided, the start of the signal will be used as the start of the trial
    `window_of_interest` 2 element array indicating the window around the stimulus index to detect spikes for
                # Units are in sample space
                # If not provided, the entire trial is used
    `fs` Sampling frequency. Default at 30000 (samples per second)
    `intra_chan_window` Duration of same spike to avoid in a single channel. Default is 6 samples (0.2ms)
    `inter_chan_window` Duration of same spike to avoid across channels. Default is 6 samples (0.2ms)

"""
function detect_spikes(
    signal_matrix::AbstractMatrix, 
    pk_bkgnds::AbstractVector, 
    extract_window = [10,10]; 
    file_lengths::AbstractVector = [size(signal_matrix, 2)],  
    stim_inds::AbstractVector = [1], 
    window_of_interest = nothing, 
    intra_chan_window=6,  
    inter_chan_window=6 
    )
    # TODO:     for a single channel: signal_matrix can be a single vector; 
    #           pk_bkgnds can be a single float

    num_channels = size(signal_matrix,1)
    sig_len = size(signal_matrix, 2)
    # onems = round(Int, fs/1000) # WOI provided in units of sample space

    largest_peaks = zeros(sig_len+1)
    cum_len = 0
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
                    signal = signal_matrix[ch_ind, cum_len+1:cum_len+fl]
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
    return pk_inds_new, concat_sigs, pk_hts #? AF: return pk_inds_new rather than pk_inds
end

function detect_spikes(
    signal_matrix::AbstractMatrix, 
    pk_bkgnds::Float64, 
    extract_window = [10,10]; 
    file_lengths::AbstractVector = [size(signal_matrix, 2)], 
    stim_inds::AbstractVector = [1], 
    window_of_interest = nothing,
    intra_chan_window=6,
    inter_chan_window=6
    )
    
    num_channels = size(signal_matrix,1)
    sig_len = size(signal_matrix, 2)
    # onems = round(Int, fs/1000) # WOI provided in units of sample space

    largest_peaks = zeros(sig_len+1)
    cum_len = 0
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
                    signal = signal_matrix[ch_ind, cum_len+1:cum_len+fl]
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



"""
    function cluster_spikes(spike_waveforms, k_override=0, apply_pca=true, PCA_Model=nothing;
        PCA_ratio = 0.9)

    Clusters spikes using a clustering method and returns:
    * Clustering Result (centers are the Templates)
    * PCA_Model

    Uses global parameters:
    * MAX_CLUSTERS (NB: these have been moved to kwards in one method)
    * CLUSTERING_EPOCS (NB: these have been moved to kwards in one method)
    * SHOW_CLUSTER_RESULT
    * COLOURS

    k_override=0 means that the best k value is calulated based on CLUSTERING_EPOCS.
    k_override>0 (Int) will force a k value. 

    PCA_Model: if you give PCA_Model any value it will cluster spikes without perfroming PCA. Leaving this term out will generate a model 

    PCA_ratio = x::float, means take the number of PCs that explain x% of the data
    PCA_scale = true scales PCs by their contribtion to variance
"""
function cluster_spikes(spike_waveforms; k_override::Int=0, apply_pca=true, PCA_Model=nothing, PCA_ratio::Float64 = 0.9, MAX_CLUSTERS = 10, CLUSTERING_EPOCS = 100, SHOW_CLUSTER_RESULT = true)

    if apply_pca == true
        if isnothing(PCA_Model)
            PCA_Model = fit(PCA, spike_waveforms; pratio=PCA_ratio)
        end
        cluster_data = MultivariateStats.transform(PCA_Model, spike_waveforms)
    else
        cluster_data = spike_waveforms
    end

    if k_override == 0
        # global MAX_CLUSTERS
        # global CLUSTERING_EPOCS
        best_kcr, best_k_count = apply_kmeans(cluster_data, max_clusters = MAX_CLUSTERS,
            num_epocs = CLUSTERING_EPOCS)
    end

    if k_override > 0
       best_kcr = reorder_kcr(kmeans(cluster_data,k_override))
       best_k_count = k_override
    end    

    # global SHOW_CLUSTER_RESULT
    global COLOURS
    if SHOW_CLUSTER_RESULT == true
        # Figure 1) Histogram of Kmeans bestK
        f1 = figure()
        PyPlot.hist(best_k_count)
        PyPlot.title("Best Clusters over $CLUSTERING_EPOCS epocs")
        # Figure 2) PCA cluster space
        try
            show_clust3D(cluster_data, best_kcr.assignments, 1,2,3, COLOURS)
            scatter3D(best_kcr.centers[1,:], best_kcr.centers[2,:],
                best_kcr.centers[3,:], color="black", linewidths = 10)
            PyPlot.title("Cluster Solution in PCA space")
            PyPlot.xlabel("Principal Component 1")
            PyPlot.ylabel("Principal Component 2")
            PyPlot.zlabel("Principal Component 3")
        catch
            @warn "Cannot show 3D Cluster Plot"
            println(PCA_Model) # PCA has outdimensions less than 3
        end
        # Figure 3) Spike waveforms with clusters together
        pyplot_waveform(spike_waveforms, best_kcr.assignments, title="Clustered raw signals and Centres")
        if apply_pca
            templates = MultivariateStats.reconstruct(PCA_Model, best_kcr.centers)
        else
            templates = best_kcr.centers
        end
        pyplot_add_centroids(templates)
        # Figure 4) Spike waveforms with clusters in separate subplots
        pyplot_waveforms_separately(spike_waveforms, best_kcr.assignments;
            centroids=templates, title="Clustered raw signals separately")
    end

    templates = reconstruct(PCA_Model, best_kcr.centers)

    return best_kcr,templates, PCA_Model, cluster_data
end


function cluster_spikes(spike_waveforms; k_override::Int=0, apply_pca=true, PCA_Model=nothing, PCA_ratio::Float64, PCA_scale::String = "off", max_clusters = 10, num_epocs = 100, SHOW_CLUSTER_RESULT = true)

    if apply_pca == true
        if isnothing(PCA_Model)
            PCA_Model = fit(PCA, spike_waveforms; pratio=PCA_ratio)
        end

        cluster_data = MultivariateStats.transform(PCA_Model, spike_waveforms)

        if PCA_scale == "on"
        PC_scaling = principalvars(PCA_Model) / tprincipalvar(PCA_Model)
        PC_scaling = PC_scaling ./  (principalvars(PCA_Model)[1] / tprincipalvar(PCA_Model))
        cluster_data = PC_scaling .* cluster_data
        end 
    else
        cluster_data = spike_waveforms
    end

    if k_override == 0
        # global MAX_CLUSTERS
        # global CLUSTERING_EPOCS
        best_kcr, best_k_count, assign_order = apply_kmeans(cluster_data, max_clusters = max_clusters, num_epocs = num_epocs)
    end

    if k_override > 0
       best_kcr, assign_order = reorder_kcr(kmeans(cluster_data,k_override))
       best_k_count = k_override
    end    

    # global SHOW_CLUSTER_RESULT
    global COLOURS
    if SHOW_CLUSTER_RESULT == true
        # Figure 1) Histogram of Kmeans bestK
        f1 = figure()
        PyPlot.hist(best_k_count)
        PyPlot.title("Best Clusters over $num_epocs epocs")
        # Figure 2) PCA cluster space
        try
            show_clust3D(cluster_data, best_kcr.assignments, 1,2,3, COLOURS)
            scatter3D(best_kcr.centers[1,:], best_kcr.centers[2,:],
                best_kcr.centers[3,:], color="black", linewidths = 10)
            PyPlot.title("Cluster Solution in PCA space")
            PyPlot.xlabel("Principal Component 1")
            PyPlot.ylabel("Principal Component 2")
            PyPlot.zlabel("Principal Component 3")
        catch
            @warn "Cannot show 3D Cluster Plot, PCA model as less than 3 components"
            println(PCA_Model) # PCA has outdimensions less than 3
        end
        # Figure 3) Spike waveforms with clusters together
        pyplot_waveform(spike_waveforms, best_kcr.assignments, title="Clustered raw signals and Centres")
        # if apply_pca
        if PCA_scale == "on"
            templates = MultivariateStats.reconstruct(PCA_Model, best_kcr.centers ./ PC_scaling)
         else
        #     templates = best_kcr.centers
            templates = MultivariateStats.reconstruct(PCA_Model, best_kcr.centers)
        end
        pyplot_add_centroids(templates)
        # Figure 4) Spike waveforms with clusters in separate subplots
        pyplot_waveforms_separately(spike_waveforms, best_kcr.assignments;
            centroids=templates, title="Clustered raw signals separately")
    end

    templates = reconstruct(PCA_Model, best_kcr.centers)
    # templates = temps[:,assign_order]

    if apply_pca==true
        return best_kcr, templates, PCA_Model, assign_order, cluster_data  
    else
        return best_kcr, cluster_data
    end
end

"""
    function keep_and_merge_templates(templates,merges)

    Keeps or merges the templates given
    `templates` Matrix where each column is a different template
    `merges` Array of arrays such as for example:
    * [[1],[2],[3],[4], [5]] (keep all of them)
    * [[1],[2,3],[5]] (merge 2nd and 3rd template, discard 4th template)
    * [[1], [2,5], [3,4]] (merge 2nd with 5th, 3rd with 4th template)
"""
function keep_and_merge_templates(templates,merges)
    output_templates = []
    for m in merges
        push!(output_templates, mean(templates[:,m], dims=2))
    end
    return hcat(output_templates...)
end

"""
    function keep_and_merge_templates(templates,merges, assignments)

    Alternative definition to arrange assignments to merged templates

"""
function keep_and_merge_templates(templates,merges, assignments)
    output_templates = []
    output_assigns = zeros(Int,length(assignments))
    for (i,m) in enumerate(merges)
        push!(output_templates, mean(templates[:,m], dims=2))
        # Try to think of a short cut, this didn't work: output_assigns[assignments .in m] = i
        indices = []
        for i0 = 1:length(assignments)
            if assignments[i0] in m
                append!(indices, i0)
            end
        end
        output_assigns[indices] .= i
        #
    end
    return hcat(output_templates...), output_assigns
end

# function compound_templates(templates)
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
        new_centers[:,original] = kcr.centers[:,new]
    end
    return KmeansResult(new_centers, new_assigns, kcr.costs, kcr.counts, kcr.wcounts, kcr.totalcost, kcr.iterations, kcr.converged), assign_order 
end

#---------------------------------------------------
# OLD FUNCTIONS
#--------------------------------------------------


"""
    apply_kmeans(values; max_clusters = 10, num_epocs = 100)

    values::Matrix of size DxN where N is number of samples
    Returns
    `kcr_best::KmeansResult` Best KCR in the mode of the best k count
    `best_k_count::Array`

    TODO Set seeds
"""
function apply_kmeans(vals; max_clusters = 10, num_epocs = 100)
    max_clusters = min(size(vals,2), max_clusters)
    dist_mat = Distances.pairwise(Euclidean(), vals, dims=2)
    best_k_count = []
    place_holder = kmeans(vals,2)
    best_kcrs = fill(place_holder, max_clusters)
    best_silos = fill(-Inf, max_clusters)
    for epoc = 1:num_epocs
        kcrs = Vector()
        for nc = 2:max_clusters #number of clusters (nc)
            kcr = kmeans(vals, nc)
            push!(kcrs, kcr)
        end
        silos = [mean(silhouettes(kcr, dist_mat)) for kcr in kcrs]
        # println(silos) #! remove after testing
        best_silo, best_k_1 = findmax(silos)
        local best_k = best_k_1 + 1 #Because `silos` start from index 1 but `nc` starts from 2
        if best_silos[best_k_1] < best_silo
            best_kcrs[best_k] = kcrs[best_k_1]
            best_silos[best_k] = best_silo
        end
        push!(best_k_count, best_k)
    end
    best_k = mode(best_k_count)
    best_kcr = best_kcrs[best_k]
    # #  4 lines to reorder best_kcr assigments into size order for consistant plotting
    # Moved to separate function
    new_kcr, assign_order = reorder_kcr(best_kcr)   
    return new_kcr, best_k_count, assign_order
    # new_kcr = reorder_kcr(best_kcr)   
    # return new_kcr, best_k_count
end



# * apply_mapdp(values; mapdp_parameters=false)
# if mapdp_parameters == true, use global parameters
# else estimate maximal marginal distribution


"""
    function get_cluster_thresholds(values, centres, assignment; cluster_pratio = 1.0, metric=Euclidean(), x_min = 10)

    `values` size(D,N) D rows of N sample columns
    `centres` size(D,K) D rows of K clusters
    `assignment` size(N) ∈ [1..K]
    `cluster_pratio` a Float64 (default 1.0), or an array size(assignment) of pratios for each assignment
    `x_min` an integer for the minimum number an assignment must have before pratio is applied (default = 10)

    Only applies a threshold if there are more than x_min sample points, otherwise uses all data points.
"""
function get_cluster_thresholds(values, centres, assignment; cluster_pratio = 1.0, metric=Euclidean(), x_min::Int = 10)
    num_clusters = maximum(assignment)
    D = size(values,1)
    cluster_thresholds = zeros(num_clusters)
    for cid = 1:num_clusters
        cluster_values = values[:, assignment .== cid]
        cluster_size = size(cluster_values,2)
        cluster_distances = Distances.pairwise(metric, reshape(centres[:,cid],D,1), cluster_values, dims=2)[:]
        if cluster_size >= x_min
            if typeof(cluster_pratio) == Float64
                cutoff = round(Int, cluster_size * cluster_pratio)
            else
            cutoff = round(Int, cluster_size * cluster_pratio[cid])
            end
        else
            cutoff = cluster_size # All of them/last index
        end
        sorted_distances = sort(cluster_distances)
        cluster_thresholds[cid] = sorted_distances[cutoff]
    end
    return cluster_thresholds
end

# function get_cluster_thresholds(values, centres, assignment; cluster_pratio::Number = 0.9, metric=Euclidean())
#         num_clusters = maximum(assignment)
#         D = size(values,1)
#         cluster_thresholds = zeros(num_clusters)
#         for cid = 1:num_clusters
#             cluster_values = values[:, assignment .== cid]
#             cluster_size = size(cluster_values,2)
#             cluster_distances = pairwise(metric, reshape(centres[:,cid],D,1), cluster_values, dims=2)[:]
#             if cluster_size >= 10
#                 cutoff = round(Int, cluster_size * cluster_pratio)
#             else
#                 cutoff = cluster_size # All of them/last index
#             end
#             sorted_distances = sort(cluster_distances)
#             cluster_thresholds[cid] = sorted_distances[cutoff]
#         end
#         return cluster_thresholds
# end
  

"""
    function assign_points_to_clusters(values, centres, thresholds=nothing)

    Assigns the value to its closest centre if it is within threshold,
    If it is within no cluster thresholds, the value remains unassigned
    `values` size(D,N)
    `centres` size(D,K)
    `thresholds=nothing` size(K,)
"""
function assign_points_to_clusters(values, centres, thresholds=nothing; metric = Euclidean())
    if isnothing(thresholds)
        thresholds = fill(Inf, size(centres,2))
    end
    num_values = size(values,2)         # N
    num_clusters = size(centres,2)      # K
    D = size(values,1)                  # D
    assignments = zeros(Int,num_values)
    # for i = 1:num_values
    #     dists = pairwise(metric, reshape(values[:,i],D,1), centres)[:]
    #     sp = sortperm(dists)
    #     sorted_dists = dists[sp]
    #     sorted_thresh = thresholds[sp]
    #     for cid = 1:num_clusters
    #         if sorted_dists[cid] <= sorted_thresh[cid]
    #             assignments[i] = cid
    #             break
    #         end
    #     end
    # end
    for i = 1:num_values
        dists = Distances.pairwise(metric, reshape(values[:,i],D,1), centres)[:]
        closest = argmin(dists)
        if dists[closest] <= thresholds[closest]
            assignments[i] = closest
        end
    end
    return assignments
end




"""
    template_merge(input_templates, keep, merge)

    * imput_templates = an array of vectors, each vector = 1 template

    * keep = vector of indices of templates to keep [1, 2, 3]

    * merge = the indices of templates to merge (average) [[5, 4], [6, 7]].
    You can either ommit the merge argument or put [] if you have nothing to merge.
    Note this is different if you want to merge two templates, which requires each merge group
        to be enclosed: [[1,2]].

    Any template not included in either the keep or merged variables is deleted. The
    order of template output is keep then merged in the order as inputed into the function.
"""
function template_merge(input_templates, keep = [], merge = [])
    output_templates=[]
    if keep != []
        kept = hcat(input_templates[keep]...)'
        push!(output_templates, kept)
    end

    if merge != []
        for i = 1: size(merge,1)
            merged = mean(hcat(input_templates[merge[i]]...)', dims = 1 )
            push!(output_templates, merged)
        end
    end
    output_templates = vcat(output_templates...)
    return output_templates
end


# Ive modified this from the preprecessing module so it works standalone


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
function assign_to_trials(stim_inds, peak_inds, iterable, pre_trial_len, post_trial_len )
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


""" Below are functions for compound templates """



"""
    compound_vectors(vect1, vect2)

    Takes two vectors and generates all possible compound vector additions by
    sliding the two vectors across each other from one end to the other and
    taking .+ for each lag instance. First vect1 is held fixed and vect2 is
    slid, then vect2 is held fixed and vect1 is slid. The concatenation of both
    vectors is not included in the output.

    Outputs: 

    `comp_vects` - an array of all compound vectors

    `temp_indexes` - an array containing the starting index of waveform 1 and 2
    in each of their respective compound vectors

"""
function compound_vectors(vect1, vect2)
    l = length(vect1)+length(vect2)-1
    comp_vects_padded = zeros(l, l)
    mtype = Vector{ Vector{Float64} }
    comp_vects = mtype(undef, 0)
    temp_indexes = []
    display(comp_vects)
    A = vcat([vect1, zeros(length(vect2)-1)]...)
    B = vcat([zeros(length(vect1)-1), vect2]...)
    for i = 1: length(vect1)
        comp_vects_padded[i,:] = A .+ circshift(B,-i+1)
        push!(comp_vects, comp_vects_padded[i,1:l-i+1])
        push!(temp_indexes, (1,length(vect1)-i+1))
    end
    for i = 1: length(vect2)-1
        comp_vects_padded[i+length(vect1),:] = circshift(B,length(vect2)) .+ circshift(A,i)
        push!(comp_vects, comp_vects_padded[i+length(vect1),1:length(vect2)+i-1])
        push!(temp_indexes, (i+1,1))
    end
    return comp_vects, temp_indexes
end

"""
    match_waveforms(wf, comp_vect)

    Given a waveform (wf) and an array of compound vectors (comp_vect), it
    returns the index of the compound vector (num) with the lowest SqEuclidean
    distance from the waveform, as well as the number of ms at which this vector
    begins relative to the waveform (index) in order to give this minimum value.

    returns num, index
"""
function match_waveforms(wf, comp_vect, showMin = "on")
    min = []
    for c in comp_vect
        A = []
        for i in 1:(length(wf)-length(c)+1)
            push!(A, (Distances.evaluate(SqEuclidean(), c, wf[i:i+length(c)-1]))/length(c))
        end
        push!(min, findmin(A))
    end
    if showMin == "on"
        display(min)
    end
    display(findmin(min))
    return  findmin(min)[1][1], findmin(min)[1][2], findmin(min)[2] # returns distance, index and num
end

function match_all_waveforms(wf, comp_vect, threshold=50)
    matches = []       # list of all compound vector matches
                        # in format (comp_vect_number, starting_index, distance)
    for (index, c) in enumerate(comp_vect)
        for i in 1:(length(wf)-length(c)+1)
            distance = Distances.evaluate(SqEuclidean(), c, wf[i:i+length(c)-1])/length(c)
            if distance <= threshold
                push!(matches, (index, i, distance))
            end
        end
    end
    return matches
end


function plot_template_sum(wf, comp, x, y, comp_index, temp_index, distance, num, trial_num=0 )
    figure()
    candidate = vcat(zeros(comp_index-1), comp, zeros(length(wf)-length(comp)-comp_index+1))
    x_padded = vcat(zeros(temp_index[1]+comp_index-2), x, zeros(length(wf)-length(x)-comp_index))
    y_padded = vcat(zeros(temp_index[2]+comp_index-2), y, zeros(length(wf)-length(y)-comp_index))
    plot(wf)
    plot(x_padded)
    plot(y_padded)
    plot(candidate)
    PyPlot.xlabel("milliseconds (ms)")
    PyPlot.ylabel("millivolts (uV)")
    PyPlot.legend(["Original Waveform", "Template 1", "Template 2",  "Sum of Templates"])
    PyPlot.title(string("Trial ", trial_num,"\nError = ", distance, "\nCompound Template ",num))
end


"""
    deconcat_templates(concat_templates, NoChannels)

    takes concatenated templates and the number of channels and breaks down the concanenated
    templates into indivudual templates in the form of multiple matrices (one matrix for each
    unconcatenated supertemplate)
"""
function deconcat_templates(concat_templates, NoChannels)
    concat_length = size(temps_adj, 2)
    iTemp_len = convert(Int,concat_length/NoChannels)
    deconcat_temps=[]
    for i = 1: NoChannels
        push!(deconcat_temps, concat_templates[:, (i*iTemp_len-iTemp_len+1):i*iTemp_len])
    end
    return deconcat_temps
end

# x and y are the original arrays (vectors)
# c is a 2D array containing all compound vectors
# rows is the number of subplots on each figure
function plot_compound_vectors(x, y, c, rows=10)
    vec_num = size(c)[1]+2                     # number of vectors to plot
    fig_num = ceil(vec_num/rows)               # total number of figures
    for f in 1:fig_num
        # make a new figure
        fig, ax = subplots(rows, 1, sharey="all", sharex = "all", figsize=(3, 6))
        fig.suptitle("Compound Vectors")
        for i in 1:rows             # plot each row for a specific figure
            cv_index = Int64(i + (f-1)*rows)
            if cv_index <= size(c)[1]
                ax[i,1].plot(c[cv_index][:])
                PyPlot.legend(["Compound Vectors"])
            else                   # plot the x and y vectors at the end
                ax[i,1].plot(x, "r")
                ax[i+1,1].plot(y, "g")
                break
            end
        end
    end
end

"""
    rmAssign() takes `X::Array` and removes the assignment `rmAs::Int` from a vector `assignments` that refer to the columns to remove from X. 

    `Y` is the returned array with columns removed
"""
function rmAssign(X,  rmAs::Int, assigns::Vector{Int64})
    Y=X[:, map(x -> x!=rmAs, assigns)]
    return Y
end

"""
    findAssign() takes `X::Array` and selects the assignment `fnd::Int` from a vector `assignments` that refer to the columns to remove from X. 

    `Y` is the returned array with columns removed
"""
function fndAssign(X,  fnd::Int, assigns::Vector{Int64})
    Y=X[:, map(x -> x==fnd, assigns)]
    return Y
end


"""
    replace_assignments() replaces `x::Int` number of `find_val::Int` with `rep_val::Int` in vector V

"""
function replace_assignments(V, x::Int, find_val::Int, rep_val::Int=0)
    z = rand(length(V))
    y = sortperm(z)
    select_inds = findall(y .< x)
    V2 = copy(V)
    j = 1
    for i = select_inds
        if i == select_inds[j]
            if V[i] == find_val 
            V2[i] = rep_val
            end
        end
        j +=1
    end
    return V2
end