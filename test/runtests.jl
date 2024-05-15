using GeneralGraphs
using Test

g = GeneralDiGraph{Int}(() -> 1, "blogs", "blogs.txt")
g = NormalUnweightedDiGraph("blogs", "blogs.txt")
g = GeneralGraph{Int8}(() -> Int8(1), "hepth", "hepth.txt")
g = NormalUnweightedGraph("hepth", "hepth.txt")

g = GeneralDiGraph{Float64}("blogs_float", "blogs_float.txt")
g = NormalWeightedDiGraph{Float64}("blogs_float", "blogs_float.txt")
g = GeneralGraph{Float64}("hepth_float", "hepth_float.txt")
g = NormalWeightedGraph{Float64}("hepth_float", "hepth_float.txt")
