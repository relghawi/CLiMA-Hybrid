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
    Random.seed!(1234)
    return Flux.Chain(
        BatchNorm(in_dim),
        Dense(in_dim => neurons, relu),
        Dense(neurons => neurons, relu),
        Dense(neurons => out_dim,sigmoid), vshape
        )
end

function LinearHybridModel(predictors, x, out_dim, neurons, b=[1.5])
    in_dim = length(predictors)
    ch = DenseNN(in_dim, out_dim, neurons)
    LinearHybridModel(ch, predictors, x, b)
end

# let's multi dispatch

function (lhm::LinearHybridModel)(df)
    x_matrix = select_predictors(df, lhm.predictors)
    α = lhm.DenseLayers(x_matrix)
    (LAIx, p_sat, p_H2O, p_atm, LA) = select_variable(df, lhm.x)
    ŷ = α #.* (p_sat - p_H2O) ./ p_atm .* LA ### F_H2O = g_lw * (p_sat-p_H2O)/p_atm * LA ## Medlyns model
    #ŷ = α
    return (; α, ŷ)
end

function (lhm::LinearHybridModel)(df, ::Val{:infer})
    α, ŷ =  lhm(df)
    return α, ŷ
end

function (lhm::LinearHybridModel)(df, infer::Symbol)
    α, ŷ = lhm(df, Val(infer))
    return α, ŷ
end


function save_predictions_to_nc(α_list, ŷ_list, filepath::String)
    num_samples = length(α_list)
    
    # Create a DataFrame with accumulated α and ŷ values
    df = DataFrame(α=zeros(num_samples), ŷ=zeros(num_samples))
    for i in 1:num_samples
        df.α[i] = α_list[i]
        df.ŷ[i] = ŷ_list[i]
    end
    
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