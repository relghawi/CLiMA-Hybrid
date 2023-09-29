using AlgebraOfGraphics
using CairoMakie
Makie.inline!(false)
include("data.jl")
names(df)
begin
    plt = data(df) * mapping(:x2, :a_syn) * visual(Scatter, markersize=5)
    with_theme(theme_ggplot2(), resolution=(900, 600)) do
        plt |> draw
    end
end
