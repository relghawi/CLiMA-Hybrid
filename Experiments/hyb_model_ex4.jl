relative_path = "../ClimaLand_examples/"

include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64.jl")

using JLD2
import ProgressMeter: BarGlyphs, next!, Progress
using Flux
using Flux.Data: DataLoader
using FluxTraining
include("../DataUtils/losses.jl")

dict = load(joinpath(relative_path, "debug.jld2"))
wd_file =joinpath(relative_path, "clima_ds_full.nc")
df = prepare_wd(dict, wd_file)

# Specify predictors and target based on DataFrame columns
predictors = [:F_H2O]
target = :F_H2O
x = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables

trainloader, valloader, testloader, trainall, d_train, d_vali, d_test = split_data(df, target, predictors, x, batchsize=100, f=0.8, v=0.1, shuffle=true, partial=true)


# Create an instance of the LinearHybridModel
model = LinearHybridModel(predictors, x, 1, 64)

# Define the loss function using the loss function from losses.jl
loss(data) = loss(model, data)

# Train the model
opt = ADAM(0.00001,(0.9, 0.999))

function custom_loss(model, loader)
    loss = 0.0  # Initialize the loss
    for data in loader
        df, y = data
        α, ŷ = model(df, :infer)
        loss = Flux.mse(ŷ, y)
    end
    return loss
end

loss_fn = Flux.Losses.mse
opt = Flux.ADAM();
# Train the model using the custom loss function
opt = ADAM(0.00001, (0.9, 0.999))
learner = Learner(model,loss_fn)

FluxTraining.fit!(learner, 10, (trainloader, valloader))

learner = Learner(model, loss)
FluxTraining.fit!(learner, 10, (trainloader, valloader))

parameters = Flux.params(model)
for data in trainloader         
    grads = Flux.gradient(parameters) do
        df, y = data  # Extract df from the data tuple
        α, ŷ = model(df)
        loss_y = Flux.mse(ŷ, y)
        learner = Learner(model, loss_y)
        
    end   
    fit!(learner, 10, (trainloader)) 
end
test!(learner, testloader)

model_state = Flux.state(model);
jldsave(joinpath(@__DIR__,relative_path,"hybrid_clima.jld2"); model_state)



# Flux.testmode!(loss,true)

d_vali2= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
column_names = [:T_AIR, :RAD, :SWC_1, :LAIx, :p_sat, :p_H2O, :p_atm, :LA]

# Reshape the data into a 1xN matrix, where N is the number of columns
data_matrix = reshape(d_vali2, 1, :)
data_df = DataFrame(data_matrix, Symbol.(column_names))

# Save predictions to an nc file
α, ŷ = model(data_df, Val(:infer)) # bug?
save_predictions_to_nc(α, ŷ, joinpath(@__DIR__,"predictions2.nc"))