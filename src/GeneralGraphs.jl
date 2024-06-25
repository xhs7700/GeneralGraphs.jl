module GeneralGraphs

using ProgressBars, ProgressMeter
using Base.Filesystem
using StatsBase
using DataStructures
using LinearAlgebra
using SparseArrays

export AbstractGraph, GeneralGraph, NormalGraph, NormalWeightedGraph, NormalUnweightedGraph
export NormalDiGraph, GeneralDiGraph, NormalWeightedDiGraph, NormalUnweightedDiGraph
export num_nodes, num_edges
export add_node!, add_edge!
export diagadj, gremban_expanse
export LCC, LCC!
export LSCC, LSCC!

abstract type AbstractGraph end
abstract type NormalGraph <: AbstractGraph end
abstract type NormalDiGraph <: AbstractGraph end

Base.size(g::AbstractGraph) = (num_nodes(g), num_edges(g))
Base.write(file_path::AbstractString, g::AbstractGraph) = open(io -> write(io, g), file_path, "w")

include("GeneralGraph.jl")
include("GeneralDiGraph.jl")
include("NormalGraph.jl")
include("NormalDiGraph.jl")

end # module GeneralGraphs