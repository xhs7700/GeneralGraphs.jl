using GeneralGraphs
using Test

g = GeneralGraph{Int}(() -> 1, "blogs", "blogs.txt")
g = NormalUnweightedGraph("blogs", "blogs.txt")
g = GeneralDiGraph{Int8}(() -> Int8(1), "hepth", "hepth.txt")
g = NormalUnweightedDiGraph("hepth", "hepth.txt")