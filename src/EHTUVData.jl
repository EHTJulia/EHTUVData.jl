module EHTUVData
using Base
using CondaPkg
using Dates
using DataFrames
using EHTDimensionalData
using EHTUtils
using FLoops
using Formatting
using OrderedCollections
using PythonCall
using Statistics
using Unitful, UnitfulAngles, UnitfulAstro # for Units

# preserve variables for Python modules
#const numpy = PyCall.PyNULL()
#const pyfits = PyCall.PyNULL()

# Abstract UV Dataset
include("abstractuvdataset/abstract.jl")

# UVDataSet
include("uvdataset/abstract.jl")
include("uvdataset/io/uvfitsutils.jl")
include("uvdataset/io/uvfitsloader.jl")
include("uvdataset/io/uvfitswriter.jl")
include("uvdataset/scan.jl")

end
