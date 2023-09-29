using Flux
using ProgressMeter
using Random

"""
machine(df, predictors, x, target, out_dim, neurons, b;
    n_epochs=10, opt = Adam(), rseed=321, FluxBoard=nothing)

Example
    - machine(df, [:x2, :x3], :x1, :obs, 1, 5, [1.0f0]; n_epochs=10, opt = Adam(), rseed=321)

"""
function machine(df, predictors, x, target, out_dim, neurons, b;
    nepochs=500, opt = Adam(), rseed=321, FluxBoard=nothing)

    trainloader, valloader, trainall = split_data(df, target, predictors, x)
    Random.seed!(rseed)
    HybridMod = LinearHybridModel(predictors, x, out_dim, neurons, b)
    opt_state = Flux.setup(opt, HybridMod)

    p = Progress(nepochs, barglyphs=BarGlyphs("[=> ]"), color = :yellow)
    for epoch in 1:nepochs
        for d in trainloader
            ∂L∂m = gradient(loss, HybridMod, d)[1]
            Flux.update!(opt_state, HybridMod, ∂L∂m)
        end

        val_loss = loss(HybridMod, valloader |> first)
        train_loss = loss(HybridMod, trainall |> first)
        next!(p; showvalues = [
                (:epoch, epoch), (:Model_seed, rseed),
                (:neurons, neurons),
                (:validation, val_loss),
                (:training, train_loss)
                ]
            )
        if !isnothing(FluxBoard)
            l_t, l_v = FluxBoard
            push!(l_t[], Point2f(epoch, train_loss))
            push!(l_v[], Point2f(epoch, val_loss))
            notify(l_v)
            notify(l_t)
            display(current_figure(); title="HybridMod (with Makie)")
        end
    end
end