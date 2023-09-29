function loss(lhm, data)
    df, y = data
    ŷ = lhm(df, :infer)::typeof(y)
    return Flux.mse(ŷ, y)
end


function loss(lhm, data)
    df, y = data
    α, ŷ = lhm(df, :infer)
    
    # Compute the loss for both ŷ and α
    loss_y = Flux.mse(ŷ, y)
    loss_α = Flux.mse(α, zeros(length(α)))  # Adjust the target as needed
    
    # You can combine or weight the two loss terms as needed
    total_loss = loss_y #+ loss_α
    
    return total_loss
end

#data = (; df = df , y = df.obs)
#loss(HybridMod, data)
#@code_warntype loss(HybridMod, data)