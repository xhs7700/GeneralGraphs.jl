struct GeneralGraph{T<:Real} <: AbstractGraph
    name::AbstractString
    nodes::Set{Int}
    edges::Dict{Tuple{Int,Int},T}
    GeneralGraph{T}() where {T<:Real} = new{T}("", Set{Int}(), Dict{Tuple{Int,Int},T}())
    GeneralGraph{T}(name::AbstractString) where {T<:Real} = new{T}(name, Set{Int}(), Dict{Tuple{Int,Int},T}())
    GeneralGraph(name::AbstractString, nodes::Set{Int}, edges::Dict{Tuple{Int,Int},T}) where {T<:Real} = new{T}(name, nodes, edges)
    function GeneralGraph{T}(name::AbstractString, io::IO) where {T<:Real}
        println("Reading graph $name...")
        g = GeneralGraph{T}(name)
        seekend(io)
        pm = Progress(position(io); dt=0.5)
        seekstart(io)
        while !eof(io)
            line = readline(io)
            ProgressMeter.update!(pm, position(io))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            words = split(line, ('\t', ' ', ','))
            u = parse(Int, words[1])
            v = parse(Int, words[2])
            w = parse(T, words[3])
            add_edge!(g, (u, v, w))
        end
        finish!(pm)
        g
    end
    function GeneralGraph(name::AbstractString, io::IO)
        println("Reading graph $name...")
        g = GeneralGraph{Int}(name)
        seekend(io)
        pm = Progress(position(io); dt=0.5)
        seekstart(io)
        while !eof(io)
            line = readline(io)
            ProgressMeter.update!(pm, position(io))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v, w = map(x -> parse(Int, x), split(line, ('\t', ' ', ',')))
            add_edge!(g, (u, v, w))
        end
        finish!(pm)
        g
    end
    function GeneralGraph{T}(name::AbstractString, source::AbstractString) where {T<:Real}
        io = open(source, "r")
        g = GeneralGraph{T}(name, io)
        close(io)
        g
    end
    function GeneralGraph(name::AbstractString, source::AbstractString)
        io = open(source, "r")
        g = GeneralGraph(name, io)
        close(io)
        g
    end
    function GeneralGraph{T}(generator::Function, name::AbstractString, io::IO) where {T<:Real}
        println("Reading graph $name...")
        g = GeneralGraph{T}(name)
        seekend(io)
        pm = Progress(position(io); dt=0.5)
        seekstart(io)
        while !eof(io)
            line = readline(io)
            ProgressMeter.update!(pm, position(io))
            if line[begin] in "#%"
                continue
            end
            line = strip(line, ('\t', ' '))
            u, v = map(x -> parse(Int, x), split(line, ('\t', ' ', ',')))
            add_edge!(g, (u, v, generator()))
        end
        finish!(pm)
        g
    end
    function GeneralGraph{T}(generator::Function, name::AbstractString, source::AbstractString) where {T<:Real}
        io = open(source, "r")
        g = GeneralGraph{T}(generator, name, io)
        close(io)
        g
    end
end

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
    pm = Progress(length(g.edges); dt=0.5)
    for (u, v) in keys(g.edges)
        next!(pm)
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
    finish!(pm)
    root = argmin(fa)
    new_nodes = filter(x -> find(x) == root, g.nodes)
    new_edges = filter(e -> find(e.first[1]) == root && find(e.first[2]) == root, g.edges)
    return GeneralGraph("$(g.name)_LCC", new_nodes, new_edges)
end

function LCC!(g::GeneralGraph)
    println("Finding LCC on graph $(g.name)...")
    fa = DefaultDict{Int,Int}(-1)
    find(x) = fa[x] < 0 ? x : (fa[x] = find(fa[x]); fa[x])
    pm = Progress(length(g.edges); dt=0.5)
    for (u, v) in keys(g.edges)
        next!(pm)
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
    finish!(pm)
    root = argmin(fa)
    filter!(x -> find(x) == root, g.nodes)
    filter!(e -> find(e.first[1]) == root && find(e.first[2]) == root, g.edges)
end

function Base.write(io::IO, g::GeneralGraph)
    n, m = size(g)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $n Edges: $m\n")
    pm = Progress(length(g.edges); dt=0.5)
    for (u, v) in g.edges |> keys |> collect |> sort
        write(io, "$u\t$v\t$(g.edges[(u,v)])\n")
        next!(pm)
    end
    finish!(pm)
end
