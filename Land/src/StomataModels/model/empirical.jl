###############################################################################
#
# Calculate empirical gsw from the equations
#
###############################################################################
"""
    stomatal_conductance(
                model::EmpiricalStomatalModel{FT},
                leaf::Leaf{FT},
                envir::AirLayer{FT},
                β::FT
    ) where {FT<:AbstractFloat}
    stomatal_conductance(
                model::EmpiricalStomatalModel{FT},
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                β::FT
    ) where {FT<:AbstractFloat}
    stomatal_conductance(
                model::EmpiricalStomatalModel{FT},
                canopyi::CanopyLayer{FT},
                envir::AirLayer{FT},
                β::FT,
                ind::Int
    ) where {FT<:AbstractFloat}

Steady state gsw from empirical approach given
- `model` [`EmpiricalStomatalModel`](@ref) type empirical model parameter set
- `leaf` [`Leaf`] type struct
- `canopyi` [`CanopyLayer`](@ref) type struct
- `envir` [`AirLayer`] type struct
- `β` Correction factor over the g1 part of an empirical model
- `ind` Nth leaf in the canopy layer
"""
function stomatal_conductance(model::ESMBallBerry{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = leaf;
    @unpack p_atm, RH = envir;

    return g0 + g1 * RH * p_atm * FT(1e-6) * β * An / p_s
end




function stomatal_conductance(model::ESMBallBerry{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = canopyi;
    @unpack p_atm, RH = envir;

    return g0 .+ g1 * RH * p_atm * FT(1e-6) * β .* An ./ p_s
end




function stomatal_conductance(model::ESMBallBerry{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1    = model;
    @unpack An, p_s   = canopyi;
    @unpack p_atm, RH = envir;

    return g0 + g1 * RH * p_atm * FT(1e-6) * β * An[ind] / p_s[ind]
end




function stomatal_conductance(model::ESMGentine{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = leaf;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An / p_i
end




function stomatal_conductance(model::ESMGentine{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = canopyi;
    @unpack p_atm   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) * β .* An ./ p_i
end




function stomatal_conductance(model::ESMGentine{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1  = model;
    @unpack An, p_i = canopyi;
    @unpack p_atm   = envir;

    return g0 + g1 * p_atm * FT(1e-6) * β * An[ind] / p_i[ind]
end




function stomatal_conductance(model::ESMLeuning{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack d0, g0, g1             = model;
    @unpack An, p_s, p_sat, Γ_star = leaf;
    @unpack p_atm, p_H₂O           = envir;

    return g0 + g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) * β * An / (p_s - Γ_star)
end




function stomatal_conductance(model::ESMLeuning{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack d0, g0, g1     = model;
    @unpack An, p_s, p_sat = canopyi;
    @unpack Γ_star         = canopyi.ps;
    @unpack p_atm, p_H₂O   = envir;

    return g0 .+ g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) * β .* An ./ (p_s .- Γ_star)
end




function stomatal_conductance(model::ESMLeuning{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack d0, g0, g1     = model;
    @unpack An, p_s, p_sat = canopyi;
    @unpack Γ_star         = canopyi.ps;
    @unpack p_atm, p_H₂O   = envir;

    return g0 + g1 * p_atm * FT(1e-6) / (1 + (p_sat - p_H₂O)/d0) * β * An[ind] / (p_s[ind] - Γ_star)
end




function stomatal_conductance(model::ESMMedlyn{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = leaf;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β * An * FT(1.6)
end




function stomatal_conductance(model::ESMMedlyn{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 .+ p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β .* An * FT(1.6)
end




function stomatal_conductance(model::ESMMedlyn{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β * An[ind] * FT(1.6)
end

relative_path = "../../../../Hybrid_model/src/"

include(joinpath(@__DIR__, relative_path, "prep_data.jl"))
include(joinpath(@__DIR__, relative_path, "data.jl"))
include(joinpath(@__DIR__, relative_path, "model_64.jl"))
include(joinpath(@__DIR__, relative_path, "losses.jl"))


using Flux, JLD2

const predictors = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]
const x = [:LAIx, :p_sat, :p_H2O, :p_atm, :LA] # Assuming as independent variables
const hybrid_model = LinearHybridModel(predictors, x, 1, 64)

const model_state_path = joinpath(@__DIR__, relative_path, "hybrid_clima.jld2")
const model_state = JLD2.load(model_state_path, "model_state")
Flux.loadmodel!(hybrid_model, model_state)


function stomatal_conductance(model::ESMMedlyn_hybrid{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    d_vali2= [p_sat, p_sat, p_sat, p_sat, p_sat, p_sat, p_sat, p_sat]
    column_names = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]

    # Reshape the data into a 1xN matrix, where N is the number of columns
    data_matrix = reshape(d_vali2, 1, :)
    data_df = DataFrame(data_matrix, Symbol.(column_names))

    α, ŷ = hybrid_model(data_df, Val(:infer))
    α =α[1]

    return α
end
