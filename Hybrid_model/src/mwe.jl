include("data.jl")
include("model.jl")
include("losses.jl")
include("machine.jl")
include("plot_live.jl")

l_t = Observable(Point2f[])
l_v = Observable(Point2f[])

fluxBoard(l_t, l_v)

machine(df, [:x2, :x3], :x1, :obs, 1, 5, [0.0f0])

trainloader, valloader, trainall = split_data(df, :obs, [:x2, :x3], :x1)

#=
HybridMod = LinearHybridModel([:x2, :x3], :x1,  1, 5, [0.0f0])
df_k = tokeyedArray(df)

data = (; df = df_k , y = df.obs)
loss(HybridMod, data)

#@code_warntype loss(HybridMod, data)
trainloader, valloader, trainall = split_data(df, :obs, [:x2, :x3], :x1)

loss(HybridMod, valloader |> first)
=#


