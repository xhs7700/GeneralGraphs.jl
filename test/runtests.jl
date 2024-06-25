using GeneralGraphs
using Test, Base.Filesystem

g = GeneralDiGraph{Int}(() -> 1, "blogs", "blogs.txt")
g = NormalUnweightedDiGraph("blogs", "blogs.txt")
g = GeneralGraph{Int8}(() -> Int8(1), "hepth", "hepth.txt")
g = NormalUnweightedGraph("hepth", "hepth.txt")

g = GeneralDiGraph{Float64}("blogs_float", "blogs_float.txt")
g = NormalWeightedDiGraph{Float64}("blogs_float", "blogs_float.txt")
g = GeneralGraph{Float64}("hepth_float", "hepth_float.txt")
g = NormalWeightedGraph{Float64}("hepth_float", "hepth_float.txt")

if !isdir("../tmp")
    mkdir("../tmp")
end
g = NormalWeightedGraph{Int}("toy4", "toy4.txt")
write("../tmp/toy4_expansion.txt", gremban_expanse(g))
g = NormalWeightedDiGraph{Int}("toy3", "toy3_raw.txt")
write("toy3.txt", g)
write("../tmp/toy3_expansion.txt", gremban_expanse(g))
