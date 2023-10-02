relative_path = "../../Land/examples/"

include("prep_data.jl")
include("data.jl")
include("model_64.jl")
include("losses.jl")
using JLD2

dict = load(joinpath(@__DIR__, relative_path, "debug.jld2"))
wd_file =joinpath(@__DIR__, relative_path, "clima_ds_full.nc")
df = prepare_wd(dict, wd_file)

# Specify predictors and target based on DataFrame columns
predictors = [:T_AIR, :RAD, :SWC_1,:LAIx, :p_sat,:p_H2O,:p_atm,:LA]
target = :F_H2O
x = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables

trainloader, valloader, trainall, d_train, d_vali = split_data(df, target, predictors, x,f = 0.8, batchsize=100, shuffle=true, partial=true)

# Create an instance of the LinearHybridModel
model = LinearHybridModel(predictors, x, 1, 64)

# Define the loss function using the loss function from losses.jl
loss(data) = loss(model, data)

# Train the model
opt = ADAM(0.00001,(0.9, 0.999))

for epoch in 1:1
    Flux.train!(loss, Flux.params(model), trainloader, opt,cb = Flux.throttle(() -> println("training"), 10)) #,cb = Flux.throttle(() -> println("training"), 10)
end

model_state = Flux.state(model);
jldsave(joinpath(@__DIR__,"hybrid_clima.jld2"); model_state)



# Flux.testmode!(loss,true)

d_vali2= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
column_names = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]

# Reshape the data into a 1xN matrix, where N is the number of columns
data_matrix = reshape(d_vali2, 1, :)
data_df = DataFrame(data_matrix, Symbol.(column_names))

# Save predictions to an nc file
α, ŷ = model(data_df, Val(:infer))
save_predictions_to_nc(α, ŷ, joinpath(@__DIR__,"predictions2.nc"))