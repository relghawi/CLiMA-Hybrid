
include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64.jl")
include("../DataUtils/losses.jl")

using Flux, JLD2

predictors = [:T_AIR,:SWC_1,:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
x = [:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables
hybrid_model = LinearHybridModel(predictors, x, 1, 128)

model_state_path = joinpath(@__DIR__, "hybrid_clima.jld2")
model_state = JLD2.load(model_state_path, "model_state")

Flux.loadmodel!(hybrid_model, model_state)

import Land.StomataModels.stomatal_conductance
using UnPack: @unpack
using Land.StomataModels: CanopyLayer, ESMMedlyn
using Land.Photosynthesis: AirLayer


# function stomatal_conductance(model::ESMMedlynHybrid{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT},swc::FT,β::FT, ind::Int) where {FT<:AbstractFloat}
#     @unpack g0, g1            = model;
#     @unpack An, p_sat, LA,LAIx,Rn  = canopyi;
#     @unpack p_a, p_atm, p_H₂O,t_air = envir;
#     vpd = max(FT(0.001), p_sat - p_H₂O);

#     d_vali2= [t_air, Rn[1], swc, LAIx[1], p_sat, p_H₂O,p_atm,LA]
#     column_names = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]

#     # Reshape the data into a 1xN matrix, where N is the number of columns
#     data_matrix = reshape(d_vali2, 1, :)
    
#     data_df = DataFrame(data_matrix, Symbol.(column_names))
    
#     α, ŷ = hybrid_model(data_df, Val(:infer))
#     α =α[1]

#     return α
# end

function stomatal_conductance(model::ESMMedlyn{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT},swc::FT, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β * An[ind] * FT(1.6)
end

import Land.SoilPlantAirContinuum.Canopy_cond_un
using UnPack: @unpack
using Land.SoilPlantAirContinuum: SPACMono
using Land.SoilPlantAirContinuum: CanopyLayer, EmpiricalStomatalModel
using Land.Photosynthesis: AirLayer


function Canopy_cond_un(spac::SPACMono{FT}, Hyb::Bool) where {FT<:AbstractFloat}
    _g_lw::FT = 0;
    _t_veg::FT = 0;
    #LAIx::FT = 0;

    if Hyb
        for _i_can in 1:spac.n_canopy
            # @show _i_can;
            _iEN = spac.envirs[_i_can];
            _iPS = spac.plant_ps[_i_can];
            #LAIx += numerical∫(_iPS.LAIx, _iPS.LAIx)
           

            vpd = max(FT(0.001), _iPS.p_sat - _iEN.p_H₂O);
            
            
            d_vali2 = [_iEN.t_air,spac.swc[1], _iPS.LAIx[1], _iPS.p_sat,_iEN.p_H₂O, _iEN.p_atm, _iPS.LA,vpd]

            column_names = [:T_AIR,:SWC_1,:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
            data_matrix = reshape(d_vali2, 1, :)
    
            data_df = DataFrame(data_matrix, Symbol.(column_names))
    
            α, ŷ = hybrid_model(data_df, Val(:infer))
            _g_lw = α[1]
            _t_veg = _g_lw .* (_iPS.p_sat - _iEN.p_H₂O) ./ _iEN.p_atm .* _iPS.LA;
    
            #_g_lw += numerical∫(_iPS.g_lw, _iPS.LAIx);
        end;
        return  _g_lw, _t_veg  #/ spac.ga
     
    end 
end



using PkgUtility: numerical∫
# using Flux, JLD2
# include("../DataUtils/prep_data.jl")
# include("../DataUtils/data.jl")
# include("../HybridModels/model_64_gc.jl")
# include("../DataUtils/losses.jl")

# predictors = [:T_AIR,:SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
# x = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:ga_spac] # Assuming as independent variables
# hybrid_model = LinearHybridModel(predictors, x, 1, 128)

# model_state_path = joinpath(@__DIR__, "hybrid_clima_gc.jld2")
# model_state = JLD2.load(model_state_path, "model_state")

# Flux.loadmodel!(hybrid_model, model_state)

function Canopy_cond(spac::SPACMono{FT}, Hyb::Bool) where {FT<:AbstractFloat}
    _g_lw::FT = 0;
    _t_veg::FT = 0;
    LAIx::FT = 0;
    Rad::FT = 0;
    if Hyb
        for _i_can in 1:spac.n_canopy
            # @show _i_can;
            _iEN = spac.envirs[_i_can];
            _iPS = spac.plant_ps[_i_can];


            LAIx += numerical∫(_iPS.LAIx, _iPS.LAIx)
            Rad += numerical∫(_iPS.Rn, _iPS.LAIx)

            vpd = max(FT(0.001), _iPS.p_sat - _iEN.p_H₂O);
            
            
            d_vali2 = [_iEN.t_air,spac.swc[1], LAIx, _iPS.p_sat,_iEN.p_H₂O, _iEN.p_atm, _iPS.LA,vpd,spac.ga]
            column_names = [:T_AIR,:SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd,:ga_spac]
            data_matrix = reshape(d_vali2, 1, :)
    
            data_df = DataFrame(data_matrix, Symbol.(column_names))
    
            α, ŷ = hybrid_model(data_df, Val(:infer))
            _g_lw = α[1]
            _t_veg = (_g_lw .* (_iPS.p_sat - _iEN.p_H₂O) ./ _iEN.p_atm .* _iPS.LA) ./ spac.ga;
    
            #_g_lw += numerical∫(_iPS.g_lw, _iPS.LAIx);
        end;
        return  _g_lw, _t_veg   #/ spac.ga
     
    end 
end


import Land.StomataModels.prognostic_gsw!
using UnPack: @unpack
using Land.StomataModels: CanopyLayer, EmpiricalStomatalModel
using Land.Photosynthesis: AirLayer

function prognostic_gsw!(clayer::CanopyLayer{FT}, envir::AirLayer{FT}, sm::EmpiricalStomatalModel{FT},  swc::FT, β::FT, Δt::FT, Hyb::Bool) where {FT<:AbstractFloat}
    # unpack values
    @unpack g_bc, g_bw, g_lc, g_lw, g_m, g_sc, g_sw, n_leaf = clayer;
    @unpack An, p_sat, LA, LAIx, Rn = clayer;  # Updated the variable to be clayer
    @unpack p_a, p_atm, p_H₂O, t_air,wind = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);


    if Hyb
        d_vali2 = [t_air, Rn[1], swc, LAIx[1], p_sat, p_H₂O, p_atm, LA,wind,vpd]
        column_names = [:T_AIR, :RAD, :SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:WIND,:vpd]
        data_matrix = reshape(d_vali2, 1, :)
        # println("LAIx prog")
        # println(LAIx)

        data_df = DataFrame(data_matrix, Symbol.(column_names))

        α, ŷ = hybrid_model(data_df, Val(:infer))
        α = α[1]

    # update g_sw
        for iLF in 1:n_leaf
            gsw_ss = max(0, stomatal_conductance(sm, clayer, envir, swc, β, iLF))
            g_sw[iLF] += (gsw_ss - g_sw[iLF]) / clayer.τ_esm * Δt

        # update g_lw, gsc, and g_lc as well
            g_lw[iLF] = α

            g_sc[iLF] = g_sw[iLF] / FT(1.6)
            g_lc[iLF] = 1 / (1 / g_sc[iLF] + 1 / g_m[iLF] + 1 / g_bc[iLF])
        end
    else

        for iLF in 1:n_leaf
            gsw_ss = max(0, stomatal_conductance(sm, clayer, envir, swc, β, iLF))
            g_sw[iLF] += (gsw_ss - g_sw[iLF]) / clayer.τ_esm * Δt

        # update g_lw, gsc, and g_lc as well
            g_lw[iLF] = 1 / (1 / g_sw[iLF] + 1 / g_bw[iLF])
            g_sc[iLF] = g_sw[iLF] / FT(1.6)
            g_lc[iLF] = 1 / (1 / g_sc[iLF] + 1 / g_m[iLF] + 1 / g_bc[iLF])
        end
    end

    return nothing
end

