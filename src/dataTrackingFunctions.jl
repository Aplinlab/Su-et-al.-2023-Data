# Version 1.0
# * Changelog

# ! Errors

#   TODO List
    # Check boundaries for defining iDC low and iDC high (line 80)

# ? Suggestions & Queries



# Packages and working directory
using CSV
using DataFrames



"""
Returns a dictionary of DataFrames containing only the subset of specified trials matching the specified cells and criteria.

cellList - list of cells (e.g. cells with evoked activity) \n
searchDF - DataFrame generated from R analysis \n
colHeading - column name in searchDF corresponding to the relevant criteria \n
searchText - value which the function will search for in the specified column \n
"""
function filterMatchingCells(cellList, searchDF, colHeading, searchText)
    # Create empty Dictionary to be populated and returned
    matchingTrialsFull = Dict{String, DataFrame}([])

    # Filter the specified list to include only trials which match the given criteria, sorted according to cell
    for cell in cellList
        matchingTrialsCell = filter(i->(i["neuronID"] == cell && i[colHeading] == searchText), searchDF)   # define DataFrame containing subset which matches
        matchingTrialsFull[cell] = matchingTrialsCell                                                   # push subset DataFrame to output dictionary with the cell name as its key
    end
    return matchingTrialsFull
end



"""
Returns a dictionary of dictionaries such that the desired counts can be accessed using dict[stim][cell].
"""
function createCountDicts(stimList, cellList, searchDF, colHeading)
    matchingCellsFull = Dict{String, Dict{String, Int64}}([])       # empty output dictionary
    for stim in stimList
        matchingCollectedStim = filterMatchingCells(cellList, searchDF, colHeading, stim)               # create dictionary of DataFrames
        countCollectedStim = Dict{String, Int64}([])                                                    # create empty dictionary
        for (cell, matchingTrials) in matchingCollectedStim
            countCollectedStim[cell] = nrow(matchingTrials)         # for each cell, store the number of matching trials into a dictionary
        end
        matchingCellsFull[stim] = countCollectedStim                # collate everything into output dictionary
    end
    return matchingCellsFull
end



"""
Returns a DataFrame of trial counts when given input DataFrames for collected, evoked, and iDC data cells.
"""
function createDataSummary(collectedDF, evokedDF, idcDF)
    # Create lists necessary for calling functions
    evokedCellList = sort(unique(evokedDF[!, "neuronID"]))          # list of all evoked cells
    stimSearchList = sort(unique(collectedDF[!, "stim_type_display"]))       # list of all defined stimulus types

    # Create dictionaries of dictionaries such that the desired counts can be accessed using dict[collected/evoked/idcLow/idcHigh][stim][cell]
    summaryDict = Dict{String, Dict{String, Dict{String, Int64}}}([])
    summaryDict["collected"] = createCountDicts(stimSearchList, evokedCellList, collectedDF, "stim_type_display")
    summaryDict["evoked"] = createCountDicts(stimSearchList, evokedCellList, evokedDF, "stim_type_display")

    # iDC data must be separated into iDC low and iDC high, so createCountDicts can't be used
    matchingIDClow = Dict{String, Dict{String, Int64}}([])
    matchingIDChigh = Dict{String, Dict{String, Int64}}([])
    for stim in stimSearchList
        matchingIDCStim = filterMatchingCells(evokedCellList, idcDF, "stim_type_display", stim)
        countIDClowStim = Dict{String, Int64}([])
        countIDChighStim = Dict{String, Int64}([])
        for (cell, matchingTrials) in matchingIDCStim
            countIDClowStim[cell] = nrow(filter(i->(0 < i["idc_amp_num"] < 1500), matchingTrials))
            countIDChighStim[cell] = nrow(filter(i->(1500 < i["idc_amp_num"] < 10000), matchingTrials))
        end
        matchingIDClow[stim] = countIDClowStim
        matchingIDChigh[stim] = countIDChighStim
    end
    summaryDict["idcLow"] = matchingIDClow
    summaryDict["idcHigh"] = matchingIDChigh



    # Create and return a DataFrame containing summary data
    summaryData = DataFrame(cells = evokedCellList)
    for stim in stimSearchList
        for condition in ["collected", "evoked", "idcLow", "idcHigh"]
            summaryCol = []
            for cell in evokedCellList
                push!(summaryCol, summaryDict[condition][stim][cell])
            end
            summaryData[!, "$stim $condition"] = summaryCol
        end
    end
    return summaryData
end