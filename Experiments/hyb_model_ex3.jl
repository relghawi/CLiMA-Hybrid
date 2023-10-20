relative_path = "../ClimaLand_examples/"

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

trainloader, valloader, testloader, trainall, d_train, d_vali, d_test = split_data(df, target, predictors, x, batchsize=100, f=0.8, v=0.1, shuffle=true, partial=true)

function getLoss(model, loader::Flux.Data.DataLoader)
    num_batches = length(loader)
    loss::Float32 = 0.0f0
    #Flux.testmode!(model)
    for data in loader
        df, y = data  # Extract df from the data tuple
        α, ŷ = model(df)
        loss_y = Flux.mse(ŷ, y)
        loss += loss_y
    end
    return loss / num_batches
    println("loss: ", loss / num_batches)
end


function train(model, trainloader::Flux.Data.DataLoader, valloader::Flux.Data.DataLoader,
        checkpoint_file::String; optimizer=Flux.ADAM(), n_epochs= 500, patience=100)

    num_epochs_no_improvement = 0
    best_epoch = 0
    best_loss = getLoss(model, valloader)
    #Flux.trainmode!(model)
    optim = Flux.setup(optimizer, model)

    bar = Progress(
        n_epochs,
        dt=1.0;
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=20,
        showspeed=true
    )
    train_losses = Float32[]
    val_losses = Float32[]

    for epoch in 1:n_epochs
        for data in trainloader            
            loss, grads = Flux.withgradient(model) do m
                df, y = data
                α, ŷ = m(df)
                loss = Flux.mse(ŷ, y)
            end
            Flux.update!(optim, model, grads[1])
        end
        train_loss = getLoss(model, trainloader)
        push!(train_losses, train_loss)
        
        val_loss = getLoss(model, valloader)
        push!(val_losses, val_loss)

        if val_loss < best_loss
            best_loss = val_loss
            num_epochs_no_improvement = 0
            best_epoch = epoch
            # Flux.state(model) |> BSON.@save checkpoint_file 
            model_state = Flux.state(model)
            jldsave(joinpath(@__DIR__,"hybrid_clima.jld2"); model_state)
            #println("Best Model Saved")
        else
            num_epochs_no_improvement += 1
            if num_epochs_no_improvement >= patience
                println("Early stopping! No improvement for $(patience) epochs.")
                break
            end
        end
        next!(bar; showvalues = [(:epoch, epoch), (:val_loss, val_loss), (:best_epoch, best_epoch)])
    end
    
    println("Training completed!")
    return train_losses, val_losses
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
checkpoint_file = "hybrid_clima.jld2"  # Specify the checkpoint file path here
train_losses, val_losses = train(model, trainloader, valloader, checkpoint_file);

load_model(model, checkpoint_file)
predictions = test_model(model, testloader)
