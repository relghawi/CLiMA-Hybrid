relative_path = "../../Land/examples/"

include("prep_data.jl")
include("data.jl")
include("model_64.jl")

using JLD2
import ProgressMeter: BarGlyphs, next!, Progress

dict = load(joinpath(@__DIR__, relative_path, "debug.jld2"))
wd_file =joinpath(@__DIR__, relative_path, "clima_ds_full.nc")
df = prepare_wd(dict, wd_file)

# Specify predictors and target based on DataFrame columns
predictors = [:F_H2O]
target = :F_H2O
x = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA] # Assuming as independent variables

trainloader, valloader, testloader, trainall, d_train, d_vali, d_test = split_data(df, target, predictors, x, batchsize=100, f=0.8, v=0.1, shuffle=true, partial=true)

function validate(model, loader::Flux.Data.DataLoader)
    num_batches = length(loader)
    loss::Float32 = 0.0

    Flux.testmode!(model)

    for data in loader
        df, y = data  # Extract df from the data tuple
        α, ŷ = model(df, :infer)
        loss_y = Flux.mse(ŷ, y)
        loss += loss_y
    end
    return loss / num_batches
    println("loss: ", loss / num_batches)
end


function train(model, trainloader::Flux.Data.DataLoader, valloader::Flux.Data.DataLoader, optimizer, checkpoint_file::String)

    num_epochs_no_improvement = 0
    best_loss = Inf
    Flux.trainmode!(model)
    parameters = Flux.params(model)
    
    bar = Progress(
        length(trainloader),
        dt=1.0;
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=20,
        showspeed=true
    )
  
    for epoch in 1:12
        for data in trainloader            
            grads = Flux.gradient(parameters) do
                df, y = data
                α, ŷ = model(df, :infer)
                loss = Flux.mse(ŷ, y)
            end
        
            Flux.Optimise.update!(optimizer, parameters, grads)
        
            next!(bar)
        end

        if epoch % 10 == 0
            loss_val = validate(model, valloader)
            println("Validation loss after epoch $epoch: $loss_val")

            if loss_val < best_loss
                best_loss = loss_val
                num_epochs_no_improvement = 0
                # Flux.state(model) |> BSON.@save checkpoint_file 
                model_state = Flux.state(model)
                jldsave(joinpath(@__DIR__,"hybrid_clima.jld2"); model_state)
                println("Best Model Saved")
            else
                num_epochs_no_improvement += 1
                if num_epochs_no_improvement >= 10
                    println("Early stopping! No improvement for 10 epochs.")
                    break
                end
            end
        end
    end
    
    println("Training completed!")
end


# Function to load the trained model from a checkpoint file
function load_model(model, checkpoint_file::String)
    model_state_path = joinpath(@__DIR__, checkpoint_file)
    model_state = JLD2.load(model_state_path, "model_state")
    Flux.loadmodel!(model, model_state)
end

# Function to test or make predictions using the loaded model
function test_model(model, testloader::Flux.Data.DataLoader)
    Flux.testmode!(model)
    all_predictions = []
    
    for data in testloader
        df, _ = data  # Assuming you don't have labels in the test data
        α, ŷ = model(df, :infer)
        save_predictions_to_nc(α, ŷ, joinpath(@__DIR__,"predictions3.nc"))
        push!(all_predictions, ŷ)
    end
    
    return hcat(all_predictions...)
end


# Create an instance of the LinearHybridModel
model = LinearHybridModel(predictors, x, 1, 256)

# Train the model
opt = ADAM(0.001,(0.9, 0.999))


checkpoint_file = "hybrid_clima.jld2"  # Specify the checkpoint file path here
train(model, trainloader, valloader, opt, checkpoint_file) 

load_model(model, checkpoint_file)
predictions = test_model(model, testloader)
