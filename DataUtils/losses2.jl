using Flux

# Adjust the loss function to work with SPAC
function loss(spac::SPACMono{FT}, dfr::DataFrame, target_variable::Symbol) where {FT<:AbstractFloat}
    # Get the predicted values for the target variable from SPAC
    ŷ = spac(dfr, :infer)[:, Symbol(target_variable)]
    
    # Extract the true target values from the DataFrame
    y_true = dfr[!, target_variable]
    
    return Flux.mse(ŷ, y_true)
end
#data = (; df = df , y = df.obs)
#loss(HybridMod, data)
#@code_warntype loss(HybridMod, data)