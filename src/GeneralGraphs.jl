module GeneralGraphs

using ProgressBars
using Base.Filesystem
using StatsBase
using DataStructures
using LinearAlgebra
using SparseArrays

export AbstractGraph, GeneralGraph, NormalGraph, NormalWeightedGraph, NormalUnweightedGraph
export NormalDiGraph, GeneralDiGraph, NormalWeightedDiGraph, NormalUnweightedDiGraph
export num_nodes, num_edges
export add_node!, add_edge!
export diagadj
export LCC, LCC!
export LSCC, LSCC!

abstract type AbstractGraph end
abstract type NormalGraph <: AbstractGraph end
abstract type NormalDiGraph <: AbstractGraph end

struct GeneralGraph{T<:Real} <: AbstractGraph
    name::AbstractString
    nodes::Set{Int}
    edges::Dict{Tuple{Int,Int},T}
    GeneralGraph{T}() where {T<:Real} = new{T}("", Set{Int}(), Dict{Tuple{Int,Int},T}())
    GeneralGraph{T}(name::AbstractString) where {T<:Real} = new{T}(name, Set{Int}(), Dict{Tuple{Int,Int},T}())
    GeneralGraph(name::AbstractString, nodes::Set{Int}, edges::Dict{Tuple{Int,Int},T}) where {T<:Real} = new{T}(name, nodes, edges)
    function GeneralGraph(name::AbstractString, source::Union{AbstractString,IO})
        println("Reading graph $name...")
        g = GeneralGraph{Int}(name)
        for line in ProgressBar(eachline(source))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v, w = map(x -> parse(Int, x), split(line, ('\t', ' ', ',')))
            add_edge!(g, (u, v, w))
        end
        g
    end
    function GeneralGraph{T}(generator::Function, name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real}
        println("Reading graph $name...")
        g = GeneralGraph{Int}(name)
        for line in ProgressBar(eachline(source))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v = map(x -> parse(Int, x), split(line, ('\t', ' ', ',')))
            add_edge!(g, (u, v, generator()))
        end
        g
    end
end

Base.size(g::AbstractGraph) = (num_nodes(g), num_edges(g))
num_nodes(g::GeneralGraph) = length(g.nodes)
num_edges(g::GeneralGraph) = length(g.edges)
prefix(::GeneralGraph) = "GeneralUndiGraph"

add_node!(g::GeneralGraph, u::Int) = push!(g.nodes, u)
function add_edge!(g::GeneralGraph{T}, (u, v, w)::Tuple{Int,Int,T}) where {T<:Real}
    if u == v
        return
    end
    if u > v
        u, v = v, u
    end
    if haskey(g.edges, (u, v))
        return
    end
    push!(g.nodes, u, v)
    g.edges[(u, v)] = w
end

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
    return GeneralGraph("$(g.name)_LCC", new_nodes, new_edges)
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

Base.write(file_path::AbstractString, g::AbstractGraph) = open(io -> write(io, g), file_path, "w")

function Base.write(io::IO, g::GeneralGraph)
    n, m = size(g)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $n Edges: $m\n")
    for (u, v) in g.edges |> keys |> collect |> sort |> ProgressBar
        write(io, "$u\t$v\t$(g.edges[(u,v)])\n")
    end
end

struct GeneralDiGraph{T<:Real} <: AbstractGraph
    name::AbstractString
    adjs::Dict{Int,Vector{Tuple{Int,T}}}
    edges::Dict{Tuple{Int,Int},T}
    GeneralDiGraph{T}() where {T<:Real} = new{T}("", Dict{Int,Vector{Tuple{Int,T}}}(), Dict{Tuple{Int,Int},T}())
    GeneralDiGraph{T}(name::AbstractString) where {T<:Real} = new{T}(name, Dict{Int,Vector{Tuple{Int,T}}}(), Dict{Tuple{Int,Int},T}())
    GeneralDiGraph{T}(name::AbstractString, adjs::Dict{Int,Vector{Tuple{Int,T}}}, edges::Dict{Tuple{Int,Int},T}) where {T<:Real} = new{T}(name, adjs, edges)
    function GeneralDiGraph(name::AbstractString, source::Union{AbstractString,IO})
        println("Reading directed graph $name...")
        g = GeneralDiGraph{Int}(name)
        for line in ProgressBar(eachline(source))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v, w = map(x -> parse(Int, x), split(line, ('\t', ' ', ',')))
            add_edge!(g, (u, v, w))
        end
        g
    end
    function GeneralDiGraph{T}(generator::Function, name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real}
        println("Reading directed graph $name...")
        g = GeneralDiGraph{Int}(name)
        for line in ProgressBar(eachline(source))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v = map(x -> parse(Int, x), split(line, ('\t', ' ', ',')))
            add_edge!(g, (u, v, generator()))
        end
        g
    end
end

num_nodes(g::GeneralDiGraph) = length(g.adjs)
num_edges(g::GeneralDiGraph) = length(g.edges)
prefix(::GeneralDiGraph) = "GeneralDiGraph"

function add_node!(g::GeneralDiGraph{T}, u::Int) where {T<:Real}
    if haskey(g.adjs, u)
        return
    end
    g.adjs[u] = Tuple{Int,T}[]
end
function add_edge!(g::GeneralDiGraph{T}, (u, v, w)::Tuple{Int,Int,T}) where {T<:Real}
    if u == v || haskey(g.edges, (u, v))
        return
    end
    add_node!(g, u), add_node!(g, v)
    push!(g.adjs[u], (v, w))
    g.edges[(u, v)] = w
end

function Base.write(io::IO, g::GeneralDiGraph)
    n, m = size(g)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $n Edges: $m\n")
    for (u, v) in g.edges |> keys |> collect |> sort |> ProgressBar
        write(io, "$u\t$v\t$(g.edges[(u,v)])\n")
    end
end


function get_lscc_info(g::GeneralDiGraph)
    println("Finding LSCC of directed graph $(g.name)...")
    n = num_nodes(g)
    dfn, low = Dict{Int,Int}(), Dict{Int,Int}()
    scc_id = Dict{Int,Int}()
    in_stk, out_of_stk, vis = Set{Int}(), Set{Int}(), Set{Int}()
    dfs_stk, stk, scc_size = Tuple{Int,Int}[], Int[], Int[]
    dfn_cnt = 0
    scc_cnt = 0
    sizehint!(dfn, n), sizehint!(low, n)
    sizehint!(vis, n), sizehint!(scc_id, n)
    sizehint!(out_of_stk, n)
    for i in keys(g.adjs)
        if i in vis
            continue
        end
        push!(dfs_stk, (0, i))
        while !isempty(dfs_stk)
            pre, u = dfs_stk[end]
            if !(u in vis)
                push!(vis, u)
                dfn_cnt += 1
                dfn[u] = low[u] = dfn_cnt
                push!(stk, u)
                push!(in_stk, u)
                for (v, _) in g.adjs[u]
                    if !(v in vis)
                        push!(dfs_stk, (u, v))
                    elseif v in in_stk
                        low[u] = min(low[u], dfn[v])
                    end
                end
            else
                pop!(dfs_stk)
                if !(u in out_of_stk)
                    push!(out_of_stk, u)
                    if pre != 0
                        low[pre] = min(low[pre], low[u])
                    end
                    if dfn[u] == low[u]
                        scc_cnt += 1
                        push!(scc_size, 0)
                        while true
                            v = pop!(stk)
                            delete!(in_stk, v)
                            scc_id[v] = scc_cnt
                            scc_size[scc_cnt] += 1
                            if v == u
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    lscc_id = argmax(scc_size)
    (lscc_id, scc_id)
end

function LSCC(g::GeneralDiGraph{T}) where {T<:Real}
    lscc_id, scc_id = get_lscc_info(g)
    new_name = "$(g.name)_LSCC"
    new_adjs = Dict{Int,Vector{Tuple{Int,T}}}()
    for u in keys(g.adjs)
        if scc_id[u] == lscc_id
            new_adjs[u] = filter(t -> scc_id[t[1]] == lscc_id, g.adjs[u])
        end
    end
    new_edges = filter(g.edges) do p
        u, v = p.first
        scc_id[u] == lscc_id && scc_id[v] == lscc_id
    end
    GeneralDiGraph{T}(new_name, new_adjs, new_edges)
end

function LSCC!(g::GeneralDiGraph)
    lscc_id, scc_id = get_lscc_info(g)
    filter!(p -> scc_id[p.first] == lscc_id, g.adjs)
    for u in keys(g.adjs)
        filter!(t -> scc_id[t[1]] == lscc_id, g.adjs[u])
    end
    filter!(g.edges) do p
        u, v = p.first
        scc_id[u] == lscc_id && scc_id[v] == lscc_id
    end
end

struct NormalWeightedGraph{T<:Real} <: NormalGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    weights::Vector{Vector{T}}
    normalized_weights::Vector{ProbabilityWeights}
    function NormalWeightedGraph(g::GeneralGraph{T}) where {T<:Real}
        n = num_nodes(g)
        o2n = DefaultDict{Int,Int}(() -> length(o2n) + 1)
        degs = zeros(Int, n)
        adjs = Vector{Vector{Int}}(undef, n)
        weights = Vector{Vector{T}}(undef, n)
        renumber = maximum(g.nodes; init=0) != n
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
        new{T}(n, length(g.edges), g.name, adjs, weights, normalized_weights)
    end
    NormalWeightedGraph(name::AbstractString, source::Union{AbstractString,IO}) = NormalWeightedGraph(GeneralGraph(name, source))
    NormalWeightedGraph{T}(generator::Function, name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real} = NormalWeightedGraph(GeneralGraph{T}(generator, name, source))
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

function diagadj(g::NormalWeightedGraph{T}) where {T<:Real}
    n, m = g.n, g.m
    A_I, A_J, A_V = Int[], Int[], T[]
    sizehint!(A_I, 2 * m), sizehint!(A_J, 2 * m), sizehint!(A_V, 2 * m)
    d = Vector{T}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i]), append!(A_V, g.weights[i])
        d[i] = sum(g.weights[i])
    end
    return d, sparse(A_I, A_J, A_V)
end

struct NormalUnweightedGraph <: NormalGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    function NormalUnweightedGraph(g::GeneralGraph)
        n = num_nodes(g)
        o2n = DefaultDict{Int,Int}(() -> length(o2n) + 1)
        degs = zeros(Int, n)
        adjs = Vector{Vector{Int}}(undef, n)
        renumber = maximum(g.nodes; init=0) != n
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
    NormalUnweightedGraph(n::Int, m::Int, name::AbstractString, adjs::Vector{Vector{Int}}) = new(n, m, name, adjs)
    NormalUnweightedGraph(name::AbstractString, source::Union{AbstractString,IO}) = NormalUnweightedGraph(GeneralGraph{Int}(() -> 1, name, source))
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

function diagadj(::Type{T}, g::NormalUnweightedGraph) where {T<:Real}
    n, m = g.n, g.m
    A_I, A_J = Int[], Int[]
    sizehint!(A_I, 2 * m), sizehint!(A_J, 2 * m)
    d = Vector{T}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i])
        d[i] = len
    end
    return d, sparse(A_I, A_J, ones(T, length(A_I)))
end

diagadj(g::NormalUnweightedGraph) = diagadj(Float64, g)

struct NormalWeightedDiGraph{T<:Real} <: NormalDiGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    weights::Vector{Vector{T}}
    normalized_weights::Vector{ProbabilityWeights}
    function NormalWeightedDiGraph(g::GeneralDiGraph{T}) where {T<:Real}
        n = num_nodes(g)
        o2n = DefaultDict{Int,Int}(() -> length(o2n) + 1)
        adjs = Vector{Vector{Int}}(undef, n)
        weights = Vector{Vector{T}}(undef, n)
        renumber = maximum(keys(g.adjs); init=0) != n
        if renumber
            for u in ProgressBar(keys(g.adjs))
                new_u = o2n[u]
                adjs[new_u] = map(t -> o2n[t[1]], g.adjs[u])
                weights[new_u] = map(t -> t[2], g.adjs[u])
            end
        else
            for u in ProgressBar(keys(g.adjs))
                adjs[u] = map(t -> t[1], g.adjs[u])
                weights[u] = map(t -> t[2], g.adjs[u])
            end
        end
        for i in (1:n)
            p = sortperm(adjs[i])
            permute!(adjs[i], p), permute!(weights[i], p)
        end
        normalized_weights = [ProbabilityWeights(normalize(weights[u], 1)) for u in (1:n)]
        new{T}(n, length(g.edges), g.name, adjs, weights, normalized_weights)
    end
    NormalWeightedDiGraph(name::AbstractString, source::Union{AbstractString,IO}) = NormalWeightedDiGraph(GeneralDiGraph(name, source))
    NormalWeightedDiGraph{T}(generator::Function, name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real} = NormalWeightedDiGraph(GeneralDiGraph{T}(generator, name, source))
end

num_nodes(g::NormalWeightedDiGraph) = g.n
num_edges(g::NormalWeightedDiGraph) = g.m
prefix(::NormalWeightedDiGraph) = "WeightedDiGraph"

function Base.write(io::IO, g::NormalWeightedDiGraph)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $(g.n) Edges: $(g.m)\n")
    for u in 1:g.n
        for (v, w) in zip(g.adjs[u], g.weights[u])
            write(io, "$u\t$v\t$w\n")
        end
    end
end

function diagadj(g::NormalWeightedDiGraph{T}) where {T<:Real}
    n, m = g.n, g.m
    A_I, A_J, A_V = Int[], Int[], T[]
    sizehint!(A_I, m), sizehint!(A_J, m), sizehint!(A_V, m)
    d = Vector{T}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i]), append!(A_V, g.weights[i])
        d[i] = sum(g.weights[i])
    end
    return d, sparse(A_I, A_J, A_V, n, n)
end

struct NormalUnweightedDiGraph <: NormalDiGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    function NormalUnweightedDiGraph(g::GeneralDiGraph)
        n = num_nodes(g)
        o2n = DefaultDict{Int,Int}(() -> length(o2n) + 1)
        adjs = Vector{Vector{Int}}(undef, n)
        renumber = maximum(keys(g.adjs); init=0) != n
        if renumber
            for u in ProgressBar(keys(g.adjs))
                new_u = o2n[u]
                adjs[new_u] = map(t -> o2n[t[1]], g.adjs[u])
                sort!(adjs[new_u])
            end
        else
            for u in ProgressBar(keys(g.adjs))
                adjs[u] = map(t -> t[1], g.adjs[u])
                sort!(adjs[u])
            end
        end
        new(n, length(g.edges), g.name, adjs)
    end
    NormalUnweightedDiGraph(n::Int, m::Int, name::AbstractString, adjs::Vector{Vector{Int}}) = new(n, m, name, adjs)
    NormalUnweightedDiGraph(name::AbstractString, source::Union{AbstractString,IO}) = NormalUnweightedDiGraph(GeneralDiGraph{Int}(() -> 1, name, source))
end

num_nodes(g::NormalUnweightedDiGraph) = g.n
num_edges(g::NormalUnweightedDiGraph) = g.m
prefix(::NormalUnweightedDiGraph) = "UnweightedDiGraph"

function Base.write(io::IO, g::NormalUnweightedDiGraph)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $(g.n) Edges: $(g.m)\n")
    for u in 1:g.n
        for v in g.adjs[u]
            write(io, "$u\t$v\n")
        end
    end
end

function diagadj(::Type{T}, g::NormalUnweightedDiGraph) where {T<:Real}
    n, m = g.n, g.m
    A_I, A_J = Int[], Int[]
    sizehint!(A_I, m), sizehint!(A_J, m)
    d = Vector{T}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i])
        d[i] = len
    end
    return d, sparse(A_I, A_J, ones(T, length(A_I)), n, n)
end

diagadj(g::NormalUnweightedDiGraph) = diagadj(Float64, g)

end # module GeneralGraphs