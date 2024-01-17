using GeneralGraphs
using Test

g = GeneralGraph{Int}(() -> 1, "blogs", "blogs.txt")
g = NormalUnweightedGraph("blogs", "blogs.txt")
g = GeneralDiGraph{Int}(() -> 1, "hepth", "hepth.txt")
g = NormalUnweightedDiGraph("hepth", "hepth.txt")