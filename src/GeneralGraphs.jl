module GeneralGraphs

export AbstractGraph, GeneralGraph, NormalGraph, NormalWeightedGraph, NormalUnweightedGraph

abstract type AbstractGraph end
abstract type NormalGraph <: AbstractGraph end

struct GeneralGraph <: AbstractGraph
end

struct NormalWeightedGraph <: NormalGraph
end

struct NormalUnweightedGraph <: NormalGraph
end

end # module GeneralGraphs
