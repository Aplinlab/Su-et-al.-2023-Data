#! removed if statement line 199

"""
    spike_sort_plotting.jl
    Placing all related plotting function in one file for easier version control

"""

using Plots
using PyPlot
using Printf

colours = ["green","blue","red", "y", "orange", "k", "pink", "gray",
"olive", "orange", "#8c564b", "#e377c2"]

"""
    pyplot_clusters

    X::DxN Matrix
"""
function pyplot_clusters(X,ClusterResult::ClusteringResult, K=-1; title="",xlegend="",ylegend="",zlegend="")
    pyplot_clusters(X, ClusterResult.assignments, K, title=title, xlegend=xlegend,
                    ylegend=ylegend, zlegend=zlegend)
end


function pyplot_clusters(X,z,K=-1; title="", xlegend="", ylegend="", zlegend="") #For legacy/core
    S = size(X)
    if length(S) == 1
        D = 1
        N = S[1]
    elseif length(S) == 2
        D = S[1]
        N = S[2]
    else
        throw(DimensionMismatch("Input data array dimensions " * string(S)))
    end

    figure()

    if K == -1
        K = length(Set(z))
    end

    for j = 1:K
        i = (z.==j)
        if D == 1
            X = Vector(X[:])
            cluster = X[i]
            PyPlot.scatter(cluster,cluster)
        elseif D == 2
            cluster = X[:,i]
            PyPlot.scatter(cluster[1,:],cluster[2,:])
        elseif D >=3
            cluster = X[:,i]
            PyPlot.scatter3D(cluster[1,:],cluster[2,:],cluster[3,:])
        end
    end

    PyPlot.title(title)
    PyPlot.xlabel(xlabel)
    PyPlot.ylabel(ylabel)
    PyPlot.zlabel(zlabel)
end

"""
waveform, centers are 2D matrices preferably in DxK format
"""
function pyplot_waveform(waveform, assignments, alpha=0; title="")
    S = size(waveform)
    #TODO later if waveform is not a matrix
    #DxN
    D = S[1]
    N = S[2]
    K = maximum(assignments)
    α = max(alpha, 0.5+1/N)

    # css = Plots.palette(:lightrainbow,max(K,2)).colors.colors #Must take at least 2 colours from a palette
    # colours = [(p.r, p.g, p.b, α) for p in css]         # RGBA values in a vector of tuples
    global colours
    figure()
    for k = 1:K
        PyPlot.plot(waveform[:,assignments.==k], color = colours[k])
    end
    PyPlot.title(title)
    PyPlot.xlabel("Time")
    PyPlot.ylabel("Amplitude")
end

function pyplot_add_centroids(centroids)
    D,K = size(centroids)
    # css = Plots.palette(:darkrainbow,max(K,2)).colors.colors
    # colours = [(p.r, p.g, p.b) for p in css]
    global colours
    for k = 1:K
        PyPlot.plot(centroids[:,k], color = colours[k], linewidth=5)
    end
end

"""
silhouette_values:: Array (K length) of Arrays
"""
function pyplot_silhouettes(silhouette_values::AbstractArray;title="",sort=true)
    K = length(silhouette_values)
    if sort
        sort!(silhouette_values,by=length,rev=true)
    end
    css = Plots.palette(:lightrainbow,max(K,2)).colors.colors
    colours = [(p.r, p.g, p.b) for p in css]
    figure()
    culx_end = 0
    for k = 1:K
        culx = culx_end + 1
        culx_end = culx + length(silhouette_values[k]) - 1
        PyPlot.bar(collect(culx:culx_end),silhouette_values[k],color=colours[k])
        PyPlot.xticks([])
    end
    PyPlot.xlabel("Clusters, K=$K")
    PyPlot.ylabel("Silhouette Values")
    PyPlot.title(title)
end

"""
    waveform_distance(sig_data, assignments, centres)

    sig_data:: DxN matrix
    assignments:: N Array, with Set(assignments)=1:K
    centres::DxK matrix
"""
function waveform_distance(sig_data, assignments, centres)
    cds = [] #cluster distances
    for k in Set(assignments)
        cd = mean(sum((sig_data[:, assignments.==k] .- centres[:,k]) .^ 2, dims=2))
        push!(cds,cd)
    end
    return mean(cds)
end

"""
    pyplot_waveform_6subplots(waveform, assignments, alpha=0; title="")

    `waveform` is a Matrix of size DxN
    `assignments` is an Array of length N
    `alpha` is the transparency, 0 being fully tranparent to 1 not transparent
    `centroids` is a Matrix of size DxK
    If centroids is provided, it is assumed that it is associated with the assignments
"""
function pyplot_waveform_6subplots(waveform, assignments, alpha=0; centroids=nothing,title="")
    D,N = size(waveform)
    K = maximum(assignments)
    α = max(alpha,0.5+1/N)
    css = Plots.palette(:lightrainbow,max(K÷6,2)).colors.colors #Must take at least 2 colours from a palette
    colours = [(p.r, p.g, p.b, α) for p in reverse(css)]         # RGBA values in a vector of tuples
    #Proper Julia performance indicates I should use polymorphism instead of these conditional flags
    if isnothing(centroids) 
        dcss = Plots.palette(:darkrainbow,max(K÷6,2)).colors.colors #Must take at least 2 colours from a palette
        dcolours = [(p.r, p.g, p.b) for p in reverse(dcss)]
    end
    figure()
    PyPlot.suptitle(title)
    for k = 1:K
        i = (k-1)÷6 + 1 # Colour counter
        p = k%6 # position tracker
        if p == 0
            p = 6
        end
        PyPlot.subplot(3,2,p)
        PyPlot.plot(waveform[:,assignments.==k], color = colours[i])
        if isnothing(centroids) 
            PyPlot.plot(centroids[:,k], color = dcolours[i])
        end
        PyPlot.xlabel("Time")
        PyPlot.ylabel("Amplitude")
    end
end

"""
    function pyplot_waveforms_separately(waveform, assignments, rows=3, cols=2;
        alpha=0.5, centroids=nothing, title="")

    Plots the waveforms by clusters into subplots controlled by rows*cols
    waveform:: DxN matrix
    assignments:: N length vector
    centroids:: DxK matrix or nothing
"""
function pyplot_waveforms_separately(waveform, assignments, rows=3, cols=1; alpha=0.5, centroids=nothing, title="")
    global colours #assume enough colours
    D,N = size(waveform)
    K = maximum(assignments)
    rows = K
    num_plots = rows * cols
    f,ax = PyPlot.subplots(rows, cols, sharex=true, sharey=true)
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

"""
    function pyplot_waveforms_assigned(waveforms, assignments1, assignments2;
            alpha=0.5, centroids=nothing, title="", label1="Ground Truth", label2="Template Assignment")

    waveform:: DxN matrix
    assignments::N
    centroids:: DxK matrix or nothing

    USAGE:
    pyplot_waveforms_assigned(concat_sigs_test, ground_truth, pk_assigns, centroids=templates, title=animalID)

"""
function pyplot_waveforms_assigned(waveforms, assignments1, assignments2;
        alpha=0.5, centroids=nothing, title="", label1="PCA template Assignment", label2="Temporal template Assignments")
    global colours
    K = max(maximum(assignments1), maximum(assignments2))
    f, ax = PyPlot.subplots(K+1, 2, sharex=true, sharey=true)
    f.suptitle(title)
    ax[1,1].set_ylabel("Unassigned")
    ax[K+1,1].set_xlabel(label1)
    ax[K+1,2].set_xlabel(label2)
    for p = 0:K
        c = (p == 0) ? "gray" : colours[p]
        plotA = waveforms[:, assignments1 .== p]
        plotB = waveforms[:, assignments2 .== p]
        if length(plotA) != 0
            ax[p+1, 1].plot(plotA, color = c, alpha=alpha)
        end
        if length(plotB) != 0
            ax[p+1, 2].plot(plotB, color = c, alpha=alpha)
        end
        if isnothing(centroids) && p > 0
            ax[p+1, 1].plot(centroids[:,p], color="black")
            ax[p+1, 2].plot(centroids[:,p], color="black")
        end
    end
end




"""
    function pyplot_rasters(trials_of_pks, trials_of_pk_assignments; tophigh=true)

    Plots a raster plots of the pks as dots. If the pk is assigned, it is coloured
    with the global colour scheme
    Arguments
    `trials_of_pks` A Vector of vectors which has the peak indices
    `trials_of_pk_assignments` A vector of vectors which has the peak assignments
    Keyword Arguments
    `axes` An optional parameter which allows plotting to existing graph
    `tophigh` Display option for the order of the trials, default=true meaning the last trial is displayed on the top
    `title` Title for the figure
    `stim_index` Index (in sample units) of when the stimulus happened #missing
    `fs` Sampling frequency in seconds #missing
"""
function pyplot_rasters(trials_of_pks, trials_of_pk_assignments; axes=nothing, tophigh=true, title="", stim_index=0, fs=nothing, visible_option=false)
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
end


"""
    function show_assigned_peaks(assignment, name; acc=NaN)
    Plotting function which shows waveform assignments
"""
function show_assigned_peaks(assignment, name; acc=NaN)
    global signals
    global peak_inds
    global channels
    global colours
    global concat_sig
    global temps_adj

    win_size = round(Int,size(temps_adj,2)/length(channels))
    num_ch = length(channels)

    f, ax = plt.subplots(1, 2)
    ax[1].plot(signals'); ax[1].legend(channels, title="Channels")
    for c = 1: maximum(assignment) #Setting base template
        ax[2].plot(temps_adj[c,:], color = colours[c])
    end
    ax[2].legend(colours, title = "Templates")
    for c = 1: num_ch #Dividing template into channels
        ax[2].text((c-1)/num_ch, 0, "Channel $(channels[c])", transform=ax[2].transAxes)
        ax[2].axvline((c-1)*win_size, color = "black")
    end
    for c = 1:maximum(assignment)
        local pks = peak_inds[assignment .== c]
        local csigs = concat_sig[assignment .== c,:]
        for p in pks
            ax[1].axvline(p+1,lw = 0.5, color = colours[c])
        end
        ax[2].plot(csigs', color = colours[c], alpha=0.3)
    end

    acc_s = @sprintf "%.2f" acc
    ax[1].set_title("Peaks Assigned by "*name*"\nAccuracy is:" * acc_s)
    ax[2].set_title("Waveforms in clusters by "*name)
    return f,ax
end

function show_clust3D(data, assignment, dim1, dim2, dim3, colours)
    fig = PyPlot.figure()

    if size(data,1) > 2
        for i = 1:size(assignment,1)
            PyPlot.scatter3D(data[dim1,i],data[dim2,i],data[dim3,i],color = colours[assignment[i]])
        end
    else
        println("Less than 3 PCs")
    end
end


"""
    Separates each assignment into a separate raster. It basically runs pyplot_rasters() for each assignment
"""
function separate_events(pks_plot, ass_plot)
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
    end

    for i = 1:num_assigns
        AOI = i-1
        pyplot_rasters(pksLocs[i], pksAss[i], title="Raster of assignments = $AOI")
    end
end