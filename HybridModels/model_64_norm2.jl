using Flux
using NetcdfIO: read_nc, save_nc!

include("../DataUtils/NormalizationModule.jl")
norm= NormalizationModule.MatrixNormalizationModule
########################################
# Model definition y = ax + b, where 
########################################
struct LinearHybridModel # lhm
    DenseLayers::Flux.Chain
    predictors::AbstractArray{Symbol}
    latents::AbstractArray{Symbol}
    x
    b
end
# Define the NormalizationModule


vshape(x) = reshape(x, :)

function DenseNN(in_dim, out_dim, neurons)
    Random.seed!(1234)
    return Flux.Chain(
        Dense(in_dim => neurons, relu),
        Dense(neurons => neurons, relu),
        Dense(neurons => neurons, relu),
        Dense(neurons => out_dim,hardtanh), vshape
        )
end

function LinearHybridModel(predictors,latents, x, out_dim, neurons, b=[1.5])
    in_dim = length(predictors)
    ch = DenseNN(in_dim, out_dim, neurons)
    LinearHybridModel(ch, predictors,latents, x, b)
end

# let's multi dispatch

function (lhm::LinearHybridModel)(df,norm_vars,stats_latents)
    println("normalized_data")
    println(norm_vars)
    α = lhm.DenseLayers(norm_vars)
    (LAIx, p_sat, p_H2O, p_atm, LA) = select_variable(df, lhm.x)
    α = norm.denormalize(stats_latents, α)
    # println("denormalized_α")
    # println(α)
    ŷ =  α .* (p_sat - p_H2O) ./ p_atm .* LA ### F_H2O = g_lw * (p_sat-p_H2O)/p_atm * LA ## Medlyns model
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


function save_predictions_to_nc(α_list, ŷ_list,y_list, filepath::String)
    num_samples = length(α_list)
    
    # Create a DataFrame with accumulated α and ŷ values
    df = DataFrame(α=zeros(num_samples), ŷ=zeros(num_samples),y=zeros(num_samples))
    for i in 1:num_samples
        df.α[i] = α_list[i]
        df.ŷ[i] = ŷ_list[i]
        df.y[i] = y_list[i]
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