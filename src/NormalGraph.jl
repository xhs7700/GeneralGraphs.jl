struct NormalWeightedGraph{T<:Real} <: NormalGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    weights::Vector{Vector{T}}
    normalized_weights::Vector{ProbabilityWeights}
    NormalWeightedGraph(n::Int, m::Int, name::AbstractString, adjs::Vector{Vector{Int}}, weights::Vector{Vector{T}}, normalized_weights::Vector{<:ProbabilityWeights}) where {T<:Real} = new{T}(n, m, name, adjs, weights, normalized_weights)
    function NormalWeightedGraph(g::GeneralGraph{T}) where {T<:Real}
        n = num_nodes(g)
        o2n = Dict{Int,Int}()
        degs = zeros(Int, n)
        adjs = Vector{Vector{Int}}(undef, n)
        weights = Vector{Vector{T}}(undef, n)
        renumber = maximum(g.nodes; init=0) != n
        pm = Progress(length(g.edges); dt=0.5)
        if renumber
            for u in g.nodes |> collect |> sort
                o2n[u] = length(o2n) + 1
            end
            for (u, v) in keys(g.edges)
                new_u, new_v = o2n[u], o2n[v]
                degs[new_u], degs[new_v] = degs[new_u] + 1, degs[new_v] + 1
                next!(pm)
            end
        else
            for (u, v) in keys(g.edges)
                degs[u], degs[v] = degs[u] + 1, degs[v] + 1
                next!(pm)
            end
        end
        finish!(pm)
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
        normalized_weights = [normalize(weights[u], 1) .|> abs |> pweights for u in (1:n)]
        new{T}(n, length(g.edges), g.name, adjs, weights, normalized_weights)
    end
    NormalWeightedGraph{T}(name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real} = NormalWeightedGraph(GeneralGraph{T}(name, source))
    NormalWeightedGraph(name::AbstractString, source::Union{AbstractString,IO}) = NormalWeightedGraph(GeneralGraph(name, source))
    NormalWeightedGraph{T}(generator::Function, name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real} = NormalWeightedGraph(GeneralGraph{T}(generator, name, source))
end

num_nodes(g::NormalWeightedGraph) = g.n
num_edges(g::NormalWeightedGraph) = g.m
prefix(::NormalWeightedGraph) = "WeightedUndiGraph"

function gremban_expanse(g::NormalWeightedGraph{T}) where {T<:Real}
    n, m = size(g)
    new_name = "$(g.name)_expansion"
    new_adjs = Vector{Vector{Int}}(undef, 2n)
    new_weights = Vector{Vector{T}}(undef, 2n)
    for u in 1:n
        new_adjs[u] = similar(g.adjs[u])
        new_adjs[u+n] = similar(g.adjs[u])
        new_weights[u] = similar(g.weights[u])
        new_weights[u+n] = similar(g.weights[u])
    end
    for u in 1:n
        for (i, (v, w)) in enumerate(zip(g.adjs[u], g.weights[u]))
            if w > 0
                new_adjs[u][i] = v
                new_adjs[u+n][i] = v + n
                new_weights[u][i] = new_weights[u+n][i] = w
            else
                new_adjs[u][i] = v + n
                new_adjs[u+n][i] = v
                new_weights[u][i] = new_weights[u+n][i] = -w
            end
        end
    end
    for u in 1:2n
        p = sortperm(new_adjs[u])
        permute!(new_adjs[u], p), permute!(new_weights[u], p)
    end
    new_normalized_weights = [normalize(new_weights[u], 1) |> pweights for u in (1:2n)]
    NormalWeightedGraph(2n, 2m, new_name, new_adjs, new_weights, new_normalized_weights)
end

function Base.write(io::IO, g::NormalWeightedGraph)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $(g.n) Edges: $(g.m)\n")
    pm = Progress(g.m; dt=0.5)
    for u in 1:g.n
        len = length(g.adjs[u])
        for i in 1:len
            v = g.adjs[u][i]
            if u > v
                continue
            end
            write(io, "$u\t$v\t$(g.weights[u][i])\n")
            next!(pm)
        end
    end
    finish!(pm)
end

function diagadj(g::NormalWeightedGraph{T}) where {T<:Real}
    n, m = g.n, g.m
    A_I, A_J, A_V = Int[], Int[], T[]
    sizehint!(A_I, 2 * m), sizehint!(A_J, 2 * m), sizehint!(A_V, 2 * m)
    d = Vector{T}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i]), append!(A_V, g.weights[i])
        d[i] = sum(abs, g.weights[i])
    end
    return d, sparse(A_I, A_J, A_V)
end

function diagadj(::Type{T1}, g::NormalWeightedGraph{T2}) where {T1<:Real,T2<:Real}
    n, m = g.n, g.m
    A_I, A_J, A_V = Int[], Int[], T1[]
    sizehint!(A_I, 2 * m), sizehint!(A_J, 2 * m), sizehint!(A_V, 2 * m)
    d = Vector{T1}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i]), append!(A_V, g.weights[i])
        d[i] = sum(abs, g.weights[i])
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
        o2n = Dict{Int,Int}()
        degs = zeros(Int, n)
        adjs = Vector{Vector{Int}}(undef, n)
        renumber = maximum(g.nodes; init=0) != n
        pm1 = Progress(length(g.edges); dt=0.5)
        if renumber
            for u in g.nodes |> collect |> sort
                o2n[u] = length(o2n) + 1
            end
            for (u, v) in keys(g.edges)
                new_u, new_v = o2n[u], o2n[v]
                degs[new_u], degs[new_v] = degs[new_u] + 1, degs[new_v] + 1
                next!(pm1)
            end
        else
            for (u, v) in keys(g.edges)
                degs[u], degs[v] = degs[u] + 1, degs[v] + 1
                next!(pm1)
            end
        end
        finish!(pm1)
        for i in (1:n)
            deg = degs[i]
            adjs[i] = Vector{Int}()
            sizehint!(adjs[i], deg)
        end
        pm2 = Progress(length(g.edges); dt=0.5)
        if renumber
            for (u, v) in keys(g.edges)
                new_u, new_v = o2n[u], o2n[v]
                push!(adjs[new_u], new_v), push!(adjs[new_v], new_u)
                next!(pm2)
            end
        else
            for (u, v) in keys(g.edges)
                push!(adjs[u], v), push!(adjs[v], u)
                next!(pm2)
            end
        end
        finish!(pm2)
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
    pm = Progress(g.m; dt=0.5)
    for u in 1:g.n
        len = length(g.adjs[u])
        for i in 1:len
            v = g.adjs[u][i]
            if u > v
                continue
            end
            write(io, "$u\t$v\n")
            next!(pm)
        end
    end
    finish!(pm)
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
