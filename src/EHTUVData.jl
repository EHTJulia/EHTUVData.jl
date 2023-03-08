module EHTUVData
using Base
using Dates
using Conda
using DataFrames
using DimensionalData
using EHTUtils
using FLoops
using Formatting
using OrderedCollections
using PyCall
using Statistics
using Unitful, UnitfulAngles, UnitfulAstro # for Units

# preserve variables for Python modules
const numpy = PyCall.PyNULL()
const pyfits = PyCall.PyNULL()

# Abstract UV Dataset
include("abstractuvdataset/abstract.jl")

# UVDataSet
include("uvdataset/abstract.jl")
include("uvdataset/utils.jl")
include("uvdataset/io/uvfitsutils.jl")
include("uvdataset/io/uvfitsloader.jl")
include("uvdataset/io/uvfitswriter.jl")

# UVDiskDataSet
end
