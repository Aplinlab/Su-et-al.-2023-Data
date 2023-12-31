# License
Data and code repository for Su et al. 2023
Copyright (C) 2023  Aplinlab

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License v3.0 as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Email: f.aplin@unsw.edu.au

# Main Pipeline
## loading
Use to load raw data into a JLD2 format.

## templateGen
Use to create and save template assignments for spikesorting.

## spikesorting
Use to sort recordings into cell assignments.

## csvExport
Extracts relevant data from JLD2 format into a CSV for analysis in R.

## analysisResults/...paper analysis v0.44
Main statistical analysis file.

## analysisResults/...paper analysis v0.44b
Statistical analysis of effects of DC polarity.

# Other Files and Folders
## spikesortingOutput
Files saved from spikesorting. Used as input in csvExport.

## spikesortingTemplates
Templates generated using templateGen. Used as input in spikesorting, alongside raw data files.

## savedFilesInfo.xlsx
Metadata regarding JLD2 files and spike-sorted units.

## analysisResults/fullResults.csv
Data used as input for R scripts (generated by csvExport).

## tactile_latency
Used for analysis of the distribution of peak latencies in tactile units.

## analysisResults/tactile_evoked_units.csv
Data from evoked tactile units (subset of fullResults.csv).

## analysisResults/tactile_peak_latencies.csv
Peak latencies of evoked tactile units (generated by tactile_latency).

## analysisResults/Significances.xlsx
Summary of significant differences.

## src
Functions used in the main pipeline.
