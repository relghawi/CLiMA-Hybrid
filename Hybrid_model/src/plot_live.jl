#using GLMakie
#Makie.inline!(false)
#GLMakie.set_window_config!(float=true)
function fluxBoard(l_train, l_valid)
    fig = Figure(resolution=(900, 600))
    ax = Axis(fig[1, 1], xlabel="epochs", ylabel="loss")
    lines!(ax, l_train; color=:grey70, label="training")
    lines!(ax, l_valid; color=:dodgerblue, label="validation")
    autolimits!(ax)
    fig
end