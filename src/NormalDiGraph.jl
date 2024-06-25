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
        pm = Progress(length(g.adjs); dt=0.5)
        if renumber
            for u in keys(g.adjs)
                new_u = o2n[u]
                adjs[new_u] = map(t -> o2n[t[1]], g.adjs[u])
                weights[new_u] = map(t -> t[2], g.adjs[u])
                next!(pm)
            end
        else
            for u in keys(g.adjs)
                adjs[u] = map(t -> t[1], g.adjs[u])
                weights[u] = map(t -> t[2], g.adjs[u])
                next!(pm)
            end
        end
        finish!(pm)
        for i in (1:n)
            p = sortperm(adjs[i])
            permute!(adjs[i], p), permute!(weights[i], p)
        end
        normalized_weights = [ProbabilityWeights(normalize(weights[u], 1)) for u in (1:n)]
        new{T}(n, length(g.edges), g.name, adjs, weights, normalized_weights)
    end
    NormalWeightedDiGraph{T}(name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real} = NormalWeightedDiGraph(GeneralDiGraph{T}(name, source))
    NormalWeightedDiGraph(name::AbstractString, source::Union{AbstractString,IO}) = NormalWeightedDiGraph(GeneralDiGraph(name, source))
    NormalWeightedDiGraph{T}(generator::Function, name::AbstractString, source::Union{AbstractString,IO}) where {T<:Real} = NormalWeightedDiGraph(GeneralDiGraph{T}(generator, name, source))
end

num_nodes(g::NormalWeightedDiGraph) = g.n
num_edges(g::NormalWeightedDiGraph) = g.m
prefix(::NormalWeightedDiGraph) = "WeightedDiGraph"

function Base.write(io::IO, g::NormalWeightedDiGraph)
    write(io, "# $(prefix(g)): $(g.name)\n# Nodes: $(g.n) Edges: $(g.m)\n")
    pm = Progress(g.m; dt=0.5)
    cnt = 0
    for u in 1:g.n
        for (v, w) in zip(g.adjs[u], g.weights[u])
            write(io, "$u\t$v\t$w\n")
        end
        cnt += length(g.adjs[u])
        ProgressMeter.update!(pm, cnt)
    end
    finish!(pm)
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
        pm = Progress(length(g.adjs); dt=0.5)
        if renumber
            for u in keys(g.adjs)
                new_u = o2n[u]
                adjs[new_u] = map(t -> o2n[t[1]], g.adjs[u])
                sort!(adjs[new_u])
                next!(pm)
            end
        else
            for u in keys(g.adjs)
                adjs[u] = map(t -> t[1], g.adjs[u])
                sort!(adjs[u])
                next!(pm)
            end
        end
        finish!(pm)
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
    pm = Progress(g.m; dt=0.5)
    cnt = 0
    for u in 1:g.n
        for v in g.adjs[u]
            write(io, "$u\t$v\n")
        end
        cnt += length(g.adjs[u])
        ProgressMeter.update!(pm, cnt)
    end
    finish!(pm)
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
