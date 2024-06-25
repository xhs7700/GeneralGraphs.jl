struct NormalWeightedDiGraph{T<:Real} <: NormalDiGraph
    n::Int
    m::Int
    name::AbstractString
    adjs::Vector{Vector{Int}}
    weights::Vector{Vector{T}}
    normalized_weights::Vector{ProbabilityWeights}
    NormalWeightedDiGraph(n::Int, m::Int, name::AbstractString, adjs::Vector{Vector{Int}}, weights::Vector{Vector{T}}, normalized_weights::Vector{<:ProbabilityWeights}) where {T<:Real} = new{T}(n, m, name, adjs, weights, normalized_weights)
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
        normalized_weights = [normalize(weights[u], 1) .|> abs |> pweights for u in (1:n)]
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
        d[i] = sum(abs, g.weights[i])
    end
    return d, sparse(A_I, A_J, A_V, n, n)
end

function diagadj(::Type{T1}, g::NormalWeightedDiGraph{T2}) where {T1<:Real,T2<:Real}
    n, m = g.n, g.m
    A_I, A_J, A_V = Int[], Int[], T1[]
    sizehint!(A_I, m), sizehint!(A_J, m), sizehint!(A_V, m)
    d = Vector{T1}(undef, n)
    for i in 1:n
        len = length(g.adjs[i])
        append!(A_I, repeat([i], len)), append!(A_J, g.adjs[i]), append!(A_V, g.weights[i])
        d[i] = sum(abs, g.weights[i])
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

function gremban_expanse(g::NormalWeightedDiGraph{T}) where {T<:Real}
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
    NormalWeightedDiGraph(2n, 2m, new_name, new_adjs, new_weights, new_normalized_weights)
end

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
