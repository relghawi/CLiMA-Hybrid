
include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64.jl")
include("../DataUtils/losses.jl")

using Flux, JLD2
### for example
# predictors = [:T_AIR,:Rad_in,:SWC_1,:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
# x = [:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables
# hybrid_model = LinearHybridModel(predictors, x, 1, 128)

# model_state_path = joinpath(@__DIR__, "hybrid_clima.jld2")
# model_state = JLD2.load(model_state_path, "model_state")

### example oz

predictors = [:T_AIR,:Rad_in, :SWC,:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
x = [:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables
hybrid_model = LinearHybridModel(predictors, x, 1, 128)

model_state_path = joinpath(@__DIR__, "hybrid_clima_oz.jld2")
model_state = JLD2.load(model_state_path, "model_state")


Flux.loadmodel!(hybrid_model, model_state)

import Land.StomataModels.stomatal_conductance
using UnPack: @unpack
using Land.StomataModels: CanopyLayer, ESMMedlyn
using Land.Photosynthesis: AirLayer



function stomatal_conductance(model::ESMMedlyn{FT}, canopyi::CanopyLayer{FT}, envir::AirLayer{FT},swc::FT, β::FT, ind::Int) where {FT<:AbstractFloat}
    @unpack g0, g1            = model;
    @unpack An, p_sat         = canopyi;
    @unpack p_a, p_atm, p_H₂O = envir;
    vpd = max(FT(0.001), p_sat - p_H₂O);

    return g0 + p_atm * FT(1e-6) / p_a * (1 + g1/sqrt(vpd)) * β * An[ind] * FT(1.6)
end


using UnPack: @unpack
using Land.SoilPlantAirContinuum: SPACMono
using Land.SoilPlantAirContinuum: CanopyLayer, EmpiricalStomatalModel
using Land.Photosynthesis: AirLayer


function Canopy_cond_un(spac::SPACMono{FT}, Hyb::Bool) where {FT<:AbstractFloat}
    _g_lw::FT = 0;
    _t_veg::FT = 0;

    if Hyb
        for _i_can in 1:spac.n_canopy
            _iEN = spac.envirs[_i_can];
            _iPS = spac.plant_ps[_i_can];
            _iRAD = spac.in_rad
            vpd = max(FT(0.001), _iPS.p_sat - _iEN.p_H₂O);
            
            d_vali2 = [_iEN.t_air,_iRAD.E_direct[1],spac.swc[1], _iPS.LAIx[1], _iPS.p_sat,_iEN.p_H₂O, _iEN.p_atm, _iPS.LA,vpd]

            #column_names = [:T_AIR,:Rad_in,:SWC_1,:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA,:vpd] ##example
            column_names = [:T_AIR,:Rad_in,:SWC,:LAIx_out_un, :p_sat,:p_H2O,:p_atm,:LA,:vpd] ##example oz site
            data_matrix = reshape(d_vali2, 1, :)
    
            data_df = DataFrame(data_matrix, Symbol.(column_names))
    
            α, ŷ = hybrid_model(data_df, Val(:infer))
            _g_lw = α[1]
            _t_veg = _g_lw .* (_iPS.p_sat - _iEN.p_H₂O) ./ _iEN.p_atm .* _iPS.LA;
        end

        return _g_lw, _t_veg  #/ spac.ga
     
    else

        for _i_can in 1:spac.n_canopy
            _iPS = spac.plant_ps[_i_can];
            _iEN = spac.envirs[_i_can];

            for iLF in 1:325
                _g_lw = _iPS.g_lw[iLF]
                _t_veg = _g_lw .* (_iPS.p_sat - _iEN.p_H₂O) ./ _iEN.p_atm .* _iPS.LA
            end
        end

        return _g_lw, _t_veg #/ spac.ga
    end
end




using PkgUtility: numerical∫
using Flux, JLD2
include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64_gc.jl")
include("../DataUtils/losses.jl")

## for example site
# predictors_gc = [:T_AIR,:Rad_in,:SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
# x_gc = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:ga_spac] # Assuming as independent variables
# hybrid_model_gc = LinearHybridModel_canopy(predictors_gc, x_gc, 1, 128)

# model_state_path_gc = joinpath(@__DIR__, "hybrid_clima_gc.jld2")
# model_state_gc = JLD2.load(model_state_path_gc, "model_state")

### for oz site

predictors_gc = [:T_AIR,:Rad_in,:SWC,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
x_gc = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:ga_spac] # Assuming as independent variables
hybrid_model_gc = LinearHybridModel_canopy(predictors_gc, x_gc, 1, 128)

model_state_path_gc = joinpath(@__DIR__, "hybrid_clima_gc_oz.jld2")
model_state_gc = JLD2.load(model_state_path_gc, "model_state")


Flux.loadmodel!(hybrid_model_gc, model_state_gc)

function Canopy_cond(spac::SPACMono{FT}, Hyb::Bool) where {FT<:AbstractFloat}
    _g_lw::FT = 0;
    _t_veg::FT = 0;
    LAIx::FT = 0;


    if Hyb
        for _i_can in 1:spac.n_canopy
            _iEN = spac.envirs[_i_can];
            _iPS = spac.plant_ps[_i_can];
            _iRAD = spac.in_rad

            LAIx += numerical∫(_iPS.LAIx, _iPS.LAIx)


            vpd = max(FT(0.001), _iPS.p_sat - _iEN.p_H₂O);
            
            d_vali2 = [_iEN.t_air,_iRAD.E_direct[1],spac.swc[1], LAIx, _iPS.p_sat,_iEN.p_H₂O, _iEN.p_atm, _iPS.LA,vpd,spac.ga]
            #column_names = [:T_AIR,:Rad_in,:SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd,:ga_spac] ##example
            column_names = [:T_AIR,:Rad_in,:SWC,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd,:ga_spac] ##oz site
            data_matrix = reshape(d_vali2, 1, :)
    
            data_df = DataFrame(data_matrix, Symbol.(column_names))
    
            α, ŷ = hybrid_model_gc(data_df, Val(:infer))
            _g_lw = α[1]

            _t_veg = (_g_lw .* (_iPS.p_sat - _iEN.p_H₂O) ./ _iEN.p_atm .* _iPS.LA) ./ spac.ga;
    
        end

        return _g_lw, _t_veg   #/ spac.ga
     
    else

        for _i_can in 1:spac.n_canopy
            _iPS = spac.plant_ps[_i_can];
            _iEN = spac.envirs[_i_can];
            _g_lw += numerical∫(_iPS.g_lw, _iPS.LAIx);
            _t_veg += numerical∫(_iPS.g_lw, _iPS.LAIx) * (_iPS.p_sat - _iEN.p_H₂O) / _iEN.p_atm * _iPS.LA
        end

        return _g_lw, _t_veg / spac.ga #/ spac.ga
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
        column_names = [:T_AIR, :RAD, :SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:WIND,:vpd] ##example
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

