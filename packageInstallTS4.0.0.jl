#!  If error Registry.toml is missing:
#!  1. Delete folder "C:\Users\[User]\.julia\registries\General"
#!  2. Run the command `Pkg.Registry.add("General")`
#!  Should now be able to run this file. Error occurs because Julia installation was interrupted.

# Install Julia packages
packageList = [
    "Clustering", 
    "Conda", 
    "CSV", 
    "DataFrames", 
    "Dates", 
    "Distances", 
    "DSP", 
    "FileIO", 
    "FindPeaks1D", 
    "Glob", 
    "GR", 
    "HypothesisTests", 
    "Interpolations", 
    "JLD2", 
    "LinearAlgebra", 
    "MultivariateStats", 
    "Plots", 
    "Printf", 
    "PyCall", 
    "PyPlot", 
    "Query", 
    "Statistics", 
    "StatsBase", 
    "XLSX", 
]
import Pkg
for pkg in packageList
    Pkg.add(pkg)
end

# Fix for Qt5Concurrent.dll bug
import Conda
ENV["JULIA_GR_PROVIDER"] = "GR"
Pkg.build("GR")
Conda.update()
Pkg.build("PyPlot")

# Install Python packages
import PyCall
PyCall.pyimport_conda("scipy.signal", "scipy")
PyCall.pyimport_conda("matplotlib.patheffects", "matplotlib")