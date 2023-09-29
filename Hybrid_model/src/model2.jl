using Flux
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

# Modify the forward function to receive inputs compatible with SPAC
function (lhm::LinearHybridModel)(x::AbstractArray, α::AbstractArray)
    α = lhm.DenseLayers(x)
    ŷ = α .* x .+ lhm.b
    return ŷ, α
end

function (lhm::LinearHybridModel)(x::AbstractArray)
    α = lhm.DenseLayers(x)
    ŷ = α .* x .+ lhm.b
    return ŷ
end

function (lhm::LinearHybridModel)(x::AbstractArray, ::Val{:infer})
    α = lhm.DenseLayers(x)
    ŷ = α .* x .+ lhm.b
    return ŷ
end

function (lhm::LinearHybridModel)(x::AbstractArray, infer::Symbol)
    return lhm(x, Val(infer))
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