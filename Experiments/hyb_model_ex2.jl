relative_path = "./ClimaLand_examples/"

include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64.jl")

using JLD2
import ProgressMeter: BarGlyphs, next!, Progress

dict = load(joinpath(relative_path, "debug.jld2"))
wd_file =joinpath(relative_path, "clima_ds_full.nc")
df = prepare_wd(dict, wd_file)

# Specify predictors and target based on DataFrame columns
predictors = [:F_H2O]
target = :F_H2O
x = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables

trainloader, valloader, trainall, d_train, d_vali = split_data(df, target, predictors, x,f = 0.8, batchsize=100, shuffle=true, partial=true)


function train(model, loader::Flux.Data.DataLoader, optimizer)
    num_batches = length(loader)
    Flux.trainmode!(model)
    parameters = Flux.params(model)
    bar = Progress(
        num_batches,
        dt=1.0;
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=20,
        showspeed=true
    )

    for data in loader
        grads = Flux.gradient(parameters) do
            df, y = data  # Extract df from the data tuple
            α, ŷ = model(df, :infer)
            loss_y = Flux.mse(ŷ, y)
        end
    
        Flux.Optimise.update!(optimizer, parameters, grads)
        next!(bar)
    end
end

function test(model, loader::Flux.Data.DataLoader)
    num_batches = length(loader)
    loss::Float32 = 0.0

    Flux.testmode!(model)

    for data in loader
        df, y = data  # Extract df from the data tuple
        α, ŷ = model(df, :infer)
        loss_y = Flux.mse(ŷ, y)
        loss += loss_y
    end

    println("loss: ", loss / num_batches)
end

# Create an instance of the LinearHybridModel
model = LinearHybridModel(predictors, x, 1, 256)

# Train the model
opt = ADAM(0.001,(0.9, 0.999))

for epoch in 1:1000
    train(model, trainloader, opt) #,cb = Flux.throttle(() -> println("training"), 10)
    test(model, valloader)
end

model_state = Flux.state(model);
jldsave(joinpath(@__DIR__,"hybrid_clima.jld2"); model_state)

# Flux.testmode!(loss,true)

d_vali2= [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
column_names = [:F_H2O,:LAIx, :p_sat, :p_H2O, :p_atm, :LA]

# Reshape the data into a 1xN matrix, where N is the number of columns
data_matrix = reshape(d_vali2, 1, :)
data_df = DataFrame(data_matrix, Symbol.(column_names))

# Save predictions to an nc file
α, ŷ = model(d_vali, Val(:infer))
save_predictions_to_nc(α, ŷ, joinpath(@__DIR__,"predictions3.nc"))