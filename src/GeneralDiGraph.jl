struct GeneralDiGraph{T<:Real} <: AbstractGraph
    name::AbstractString
    adjs::Dict{Int,Vector{Tuple{Int,T}}}
    edges::Dict{Tuple{Int,Int},T}
    GeneralDiGraph{T}() where {T<:Real} = new{T}("", Dict{Int,Vector{Tuple{Int,T}}}(), Dict{Tuple{Int,Int},T}())
    GeneralDiGraph{T}(name::AbstractString) where {T<:Real} = new{T}(name, Dict{Int,Vector{Tuple{Int,T}}}(), Dict{Tuple{Int,Int},T}())
    GeneralDiGraph{T}(name::AbstractString, adjs::Dict{Int,Vector{Tuple{Int,T}}}, edges::Dict{Tuple{Int,Int},T}) where {T<:Real} = new{T}(name, adjs, edges)
    function GeneralDiGraph{T}(name::AbstractString, io::IO) where {T<:Real}
        println("Reading directed graph $name...")
        g = GeneralDiGraph{T}(name)
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
    function GeneralDiGraph(name::AbstractString, io::IO)
        println("Reading directed graph $name...")
        g = GeneralDiGraph{Int}(name)
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
    function GeneralDiGraph{T}(name::AbstractString, source::AbstractString) where {T<:Real}
        io = open(source, "r")
        g = GeneralDiGraph{T}(name, io)
        close(io)
        g
    end
    function GeneralDiGraph(name::AbstractString, source::AbstractString)
        io = open(source, "r")
        g = GeneralDiGraph(name, io)
        close(io)
        g
    end
    function GeneralDiGraph{T}(generator::Function, name::AbstractString, io::IO) where {T<:Real}
        println("Reading directed graph $name...")
        g = GeneralDiGraph{T}(name)
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
    function GeneralDiGraph{T}(generator::Function, name::AbstractString, source::AbstractString) where {T<:Real}
        io = open(source, "r")
        g = GeneralDiGraph{T}(generator, name, io)
        close(io)
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
    pm = Progress(length(g.edges); dt=0.5)
    for (u, v) in g.edges |> keys |> collect |> sort
        write(io, "$u\t$v\t$(g.edges[(u,v)])\n")
        next!(pm)
    end
    finish!(pm)
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
