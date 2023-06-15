module GeneralGraphs

using ProgressBars
using Base.Filesystem
using StatsBase
using DataStructures
using LinearAlgebra

export AbstractGraph, GeneralGraph, NormalGraph, NormalWeightedGraph, NormalUnweightedGraph
export num_nodes, num_edges
export LCC, LCC!

abstract type AbstractGraph end
abstract type NormalGraph <: AbstractGraph end

struct GeneralGraph{T<:Real} <: AbstractGraph
    name::AbstractString
    nodes::Set{Int}
    edges::Dict{Tuple{Int,Int},T}
    GeneralGraph{T}() where {T<:Real} = new{T}("", Set{Int}(), Dict{Tuple{Int,Int},T}())
    GeneralGraph{T}(name::AbstractString, nodes::Set{Int}, edges::Dict{Tuple{Int,Int},T}) where {T<:Real} = new{T}(name, nodes, edges)
    function GeneralGraph(name::AbstractString, file_path::AbstractString)
        println("Reading graph from $file_path...")
        nodes = Set{Int}()
        edges = Dict{Tuple{Int,Int},Int}()
        for line in ProgressBar(eachline(file_path))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v, w = map(x -> parse(Int, x), split(line, ('\t', ' ')))
            if u == v
                continue
            end
            if u > v
                u, v = v, u
            end
            if haskey(edges, (u, v))
                continue
            end
            push!(nodes, u, v)
            edges[(u, v)] = w
        end
        new{Int}(name, nodes, edges)
    end
    function GeneralGraph{T}(generator::Function, name::AbstractString, file_path::AbstractString) where {T<:Real}
        println("Reading graph from $file_path...")
        nodes = Set{Int}()
        edges = Dict{Tuple{Int,Int},T}()
        for line in ProgressBar(eachline(file_path))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v = map(x -> parse(Int, x), split(line, ('\t', ' ')))
            if u == v
                continue
            end
            if u > v
                u, v = v, u
            end
            if haskey(edges, (u, v))
                continue
            end
            push!(nodes, u, v)
            edges[(u, v)] = generator()
        end
        new{T}(name, nodes, edges)
    end
end

Base.size(g::AbstractGraph) = (num_nodes(g), num_edges(g))
num_nodes(g::GeneralGraph) = length(g.nodes)
num_edges(g::GeneralGraph) = length(g.edges)
prefix(::GeneralGraph) = "GeneralUndiGraph"

function LCC(g::GeneralGraph)
    println("Finding LCC of graph $(g.name)...")
    fa = DefaultDict{Int,Int}(-1)
    find(x) = fa[x] < 0 ? x : (fa[x] = find(fa[x]); fa[x])
    for (u, v) in ProgressBar(keys(g.edges))
        u, v = find(u), find(v)
        if u == v
            continue
        end
        if -fa[u] < -fa[v]
            u, v = v, u
        end
        fa[u] += fa[v]
        fa[v] = u
    end
    root = argmin(fa)
    new_nodes = filter(x -> find(x) == root, g.nodes)
    new_edges = filter(e -> find(e.first[1]) == root && find(e.first[2]) == root, g.edges)
    return GeneralGraph{eltype(values(new_edges))}("$(g.name)_LCC", new_nodes, new_edges)
end

function LCC!(g::GeneralGraph)
    println("Finding LCC on graph $(g.name)...")
    fa = DefaultDict{Int,Int}(-1)
    find(x) = fa[x] < 0 ? x : (fa[x] = find(fa[x]); fa[x])
    for (u, v) in ProgressBar(keys(g.edges))
        u, v = find(u), find(v)
        if u == v
            continue
        end
        if -fa[u] < -fa[v]
            u, v = v, u
        end
        fa[u] += fa[v]
        fa[v] = u
    end
    root = argmin(fa)
    filter!(x -> find(x) == root, g.nodes)
    filter!(e -> find(e.first[1]) == root && find(e.first[2]) == root, g.edges)
end

function Base.write(io::IO, g::GeneralGraph)
    n, m = size(g)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $n Edges: $m\n")
    for (u, v) in g.edges |> keys |> collect |> sort |> ProgressBar
        write(io, "$u\t$v\t$(g.edges[(u,v)])\n")
    end
end

struct NormalWeightedGraph{T<:Real} <: NormalGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    weights::Vector{Vector{T}}
    normalized_weights::Vector{ProbabilityWeights}
    function NormalWeightedGraph{T}(g::GeneralGraph, renumber::Bool) where {T<:Real}
        n = num_nodes(g)
        o2n = DefaultDict{Int,Int}(() -> length(o2n) + 1)
        degs = zeros(Int, n)
        adjs = Vector{Vector{Int}}(undef, n)
        weights = Vector{Vector{T}}(undef, n)
        if renumber
            for (u, v) in ProgressBar(keys(g.edges))
                new_u, new_v = o2n[u], o2n[v]
                degs[new_u], degs[new_v] = degs[new_u] + 1, degs[new_v] + 1
            end
        else
            for (u, v) in ProgressBar(keys(g.edges))
                degs[u], degs[v] = degs[u] + 1, degs[v] + 1
            end
        end
        for i in (1:n)
            deg = degs[i]
            adjs[i] = Vector{Int}()
            sizehint!(adjs[i], deg)
            weights[i] = Vector{T}()
            sizehint!(weights[i], deg)
        end
        if renumber
            for ((u, v), w) in g.edges
                new_u, new_v = o2n[u], o2n[v]
                push!(adjs[new_u], new_v), push!(adjs[new_v], new_u)
                push!(weights[new_u], w), push!(weights[new_v], w)
            end
        else
            for ((u, v), w) in g.edges
                push!(adjs[u], v), push!(adjs[v], u)
                push!(weights[u], w), push!(weights[v], w)
            end
        end
        for i in (1:n)
            p = sortperm(adjs[i])
            permute!(adjs[i], p), permute!(weights[i], p)
        end
        normalized_weights = [ProbabilityWeights(normalize(weights[u], 1)) for u in (1:n)]
        new(n, length(g.edges), g.name, adjs, weights, normalized_weights)
    end
    NormalWeightedGraph(name::AbstractString, file_path::AbstractString; renumber::Bool) = NormalWeightedGraph{Int}(GeneralGraph(name, file_path), renumber)
    NormalWeightedGraph{T}(generator::Function, name::AbstractString, file_path::AbstractString; renumber::Bool) where {T<:Real} = NormalWeightedGraph{T}(GeneralGraph{T}(generator, name, file_path), renumber)
end

num_nodes(g::NormalWeightedGraph) = g.n
num_edges(g::NormalWeightedGraph) = g.m
prefix(::NormalWeightedGraph) = "WeightedUndiGraph"

function Base.write(io::IO, g::NormalWeightedGraph)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $(g.n) Edges: $(g.m)\n")
    for u in 1:g.n
        len = length(g.adjs[u])
        for i in 1:len
            v = g.adjs[u][i]
            if u > v
                continue
            end
            write(io, "$u\t$v\t$(g.weights[u][i])\n")
        end
    end
end

struct NormalUnweightedGraph <: NormalGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    function NormalUnweightedGraph(g::GeneralGraph, renumber::Bool)
        n = num_nodes(g)
        o2n = DefaultDict{Int,Int}(() -> length(o2n) + 1)
        degs = zeros(Int, n)
        adjs = Vector{Vector{Int}}(undef, n)
        if renumber
            for (u, v) in ProgressBar(keys(g.edges))
                new_u, new_v = o2n[u], o2n[v]
                degs[new_u], degs[new_v] = degs[new_u] + 1, degs[new_v] + 1
            end
        else
            for (u, v) in ProgressBar(keys(g.edges))
                degs[u], degs[v] = degs[u] + 1, degs[v] + 1
            end
        end
        for i in (1:n)
            deg = degs[i]
            adjs[i] = Vector{Int}()
            sizehint!(adjs[i], deg)
        end
        if renumber
            for (u, v) in ProgressBar(keys(g.edges))
                new_u, new_v = o2n[u], o2n[v]
                push!(adjs[new_u], new_v), push!(adjs[new_v], new_u)
            end
        else
            for (u, v) in ProgressBar(keys(g.edges))
                push!(adjs[u], v), push!(adjs[v], u)
            end
        end
        for i in (1:n)
            p = sortperm(adjs[i])
            permute!(adjs[i], p)
        end
        new(n, length(g.edges), g.name, adjs)
    end
    NormalUnweightedGraph(name::AbstractString, file_path::AbstractString; renumber::Bool) = NormalUnweightedGraph(GeneralGraph{Int}(() -> 1, name, file_path), renumber)
end

num_nodes(g::NormalUnweightedGraph) = g.n
num_edges(g::NormalUnweightedGraph) = g.m
prefix(::NormalUnweightedGraph) = "UnweightedUndiGraph"

function Base.write(io::IO, g::NormalUnweightedGraph)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $(g.n) Edges: $(g.m)\n")
    for u in 1:g.n
        len = length(g.adjs[u])
        for i in 1:len
            v = g.adjs[u][i]
            if u > v
                continue
            end
            write(io, "$u\t$v\n")
        end
    end
end

end # module GeneralGraphs
