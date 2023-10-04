using Flux
g = GRU(2=>1)
g(rand(Float32,2,15))
########################################
# Model definition y = ax + b, where 
########################################
struct GRUHybridModel # lhm
    GRU_DenseLayers::Flux.Chain
    predictors::AbstractArray{Symbol}
    x
    b
end
vshape(x) = reshape(x, :)

function GRU_NN(in_dim, out_dim, neurons)
    return Flux.Chain(
        BatchNorm(in_dim),
        GRU(in_dim=>neurons),
        Dense(neurons => out_dim), vshape
        )
end

function GRU_HybridModel(predictors, x, out_dim, neurons, b=[1.5f0])
    in_dim = length(predictors)
    ch = GRU_NN(in_dim, out_dim, neurons)
    GRU_HybridModel(ch, predictors, x, b)
end

# let's multi dispatch

function (lhm::GRU_HybridModel)(df)
    x_matrix = select_predictors(df, lhm.predictors)
    α = lhm.GRU_DenseLayers(x_matrix)
    x = select_variable(df, lhm.x)
    ŷ = α .* x .+ lhm.b
    return (; α, ŷ)
end

function (lhm::GRU_HybridModel)(df, ::Val{:infer})
    _, ŷ =  lhm(df)
    return ŷ
end

function (lhm::GRU_HybridModel)(df, infer::Symbol)
    return lhm(df, Val(infer))::Vector{Float32}
end

# Call @functor to allow for training the custom model
Flux.@functor GRU_HybridModel