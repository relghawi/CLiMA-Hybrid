relative_path = "./ClimaLand_examples/"

include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64.jl")
include("../DataUtils/losses.jl")


using Flux, JLD2

const predictors = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]
const x = [:LAIx, :p_sat, :p_H2O, :p_atm, :LA] # Assuming as independent variables
const hybrid_model = LinearHybridModel(predictors, x, 1, 64)

const model_state_path = joinpath(@__DIR__, relative_path, "hybrid_clima.jld2")
const model_state = JLD2.load(model_state_path, "model_state")
Flux.loadmodel!(hybrid_model, model_state)


function stomatal_conductance(model::ESMMedlynHybrid{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT}, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat, LA,LAIx  = canopyi;
    @unpack p_a, p_atm, p_H₂O,t_air = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);
    println(p_sat)
    d_vali2= [t_air, p_sat, p_sat, LAIx, p_sat, p_H₂O,p_atm,LA]
    column_names = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]

    # Reshape the data into a 1xN matrix, where N is the number of columns
    data_matrix = reshape(d_vali2, 1, :)
    data_df = DataFrame(data_matrix, Symbol.(column_names))

    α, ŷ = hybrid_model(data_df, Val(:infer))
    α =α[1]

    return α
end