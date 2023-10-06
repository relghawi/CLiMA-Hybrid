using Flux
using NetcdfIO: read_nc, save_nc!
########################################
# Model definition y = ax + b, where 
########################################
struct LinearHybridModel # lhm
    DenseLayers::Flux.Chain
    predictors::AbstractArray{Symbol}
    x
    b
end

vshape(x) = reshape(x, :)

function DenseNN(in_dim, out_dim, neurons)
    return Flux.Chain(
        BatchNorm(in_dim),
        Dense(in_dim => neurons, σ),
        Dense(neurons => out_dim), vshape
        )
end

function LinearHybridModel(predictors, x, out_dim, neurons, b=[1.5f0])
    in_dim = length(predictors)
    ch = DenseNN(in_dim, out_dim, neurons)
    LinearHybridModel(ch, predictors, x, b)
end

# let's multi dispatch

function (lhm::LinearHybridModel)(df)
    x_matrix = select_predictors(df, lhm.predictors)
    α = lhm.DenseLayers(x_matrix)
    x = select_variable(df, lhm.x)
    ŷ = α .* x .+ lhm.b
    return (; α, ŷ)
end

function (lhm::LinearHybridModel)(df, ::Val{:infer})
    _, ŷ =  lhm(df)
    return ŷ
end

function (lhm::LinearHybridModel)(df, infer::Symbol)
    return lhm(df, Val(infer))::Vector{Float32}
end

function save_predictions_to_nc(α, ŷ, filepath::String)
    df = DataFrame(α=α, ŷ=ŷ)
    save_nc!(filepath, df)
end
# Call @functor to allow for training the custom model
Flux.@functor LinearHybridModel


# Recurrent model def, overwriting the other (not good of course)
chain4afun(nInVar) = Flux.Chain(
    BatchNorm(nInVar, affine=true),
    GRU(nInVar => 5),
    Dense(5 => 1),
)

function GRU_NN(in_dim, out_dim, neurons)
    return Flux.Chain(
        BatchNorm(in_dim),
        GRU(in_dim=>neurons),
        Dense(neurons => out_dim), vshape
        )
end