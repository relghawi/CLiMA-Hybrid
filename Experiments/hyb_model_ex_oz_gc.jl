relative_path = "../ClimaLand_examples/"

include("../DataUtils/prep_data.jl")
include("../DataUtils/data.jl")
include("../HybridModels/model_64_gc.jl")

using JLD2
import ProgressMeter: BarGlyphs, next!, Progress

wd_file ="/Net/Groups/BGI/people/relghawi/Julia_hyb/Clima-Hyb2/ClimaLand_examples/oz_full.nc"

df= read_nc(wd_file)
df.T_AIR .= df.T_AIR .+ 273.15

df.SWC ./= 100

# Specify predictors and target based on DataFrame columns
predictors = [:T_AIR,:Rad_in,:SWC,:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:vpd]
latents=:g_lw
target = :F_H2O
x = [:LAIx, :p_sat,:p_H2O,:p_atm,:LA,:ga_spac] # Assuming as independent variables

trainloader, valloader, testloader, trainall, d_train, d_vali, d_test = split_data(df, target, predictors,latents, x, batchsize=50, f=0.8, v=0.1, shuffle=true, partial=true)


function train(model, trainloader::Flux.Data.DataLoader, valloader::Flux.Data.DataLoader,
        checkpoint_file::String; optimizer=Flux.ADAM(0.001), n_epochs=5000, patience=100)

    num_epochs_no_improvement = 0
    best_epoch = 0
    best_loss = Float64(Inf)
    #Flux.trainmode!(model)
    optim = Flux.setup(optimizer, model)

    bar = Progress(
        n_epochs,
        dt=1.0;
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=20,
        showspeed=true
    )



    for epoch in 1:n_epochs

        num_batches_train = length(trainloader)
        train_loss::Float32 = 0.0
        for data in trainloader            
            loss, grads = Flux.withgradient(model) do m
                df, y,lat = data
                α, ŷ = m(df)
                loss = Flux.mse(ŷ, y)
                train_loss += loss

            end
            Flux.update!(optim, model, grads[1])
        end
        
        num_batches_val = length(valloader)
        val_loss::Float32 = 0.0
        for data in valloader
            loss, grads = Flux.withgradient(model) do m
                df, y,lat = data       
                α, ŷ = m(df)
                loss = Flux.mse(ŷ, y)
                val_loss += loss
            end
            #println("val_loss: ", val_loss /num_batches_val )
        end
 

    
        if val_loss < best_loss
            best_loss = val_loss
            num_epochs_no_improvement = 0
            best_epoch = epoch
            # Flux.state(model) |> BSON.@save checkpoint_file 
            model_state = Flux.state(model)
            jldsave(joinpath(@__DIR__,relative_path,checkpoint_file); model_state)
            #println("Best Model Saved")
        else
            num_epochs_no_improvement += 1
            if num_epochs_no_improvement >= patience
                println("Early stopping! No improvement for $(patience) epochs.")
                break
            end
        end
        next!(bar; showvalues = [(:epoch, epoch), (:train_loss, train_loss),(:val_loss, val_loss), (:best_epoch, best_epoch)])
    end
    

    println("Training completed!")

end

# Function to load the trained model from a checkpoint file
function load_model(model, checkpoint_file::String)
    model_state_path = joinpath(@__DIR__,relative_path, checkpoint_file)
    model_state = JLD2.load(model_state_path, "model_state")
    Flux.loadmodel!(model, model_state)
end

# Function to test or make predictions using the loaded model
function test_model(model, testloader::Flux.Data.DataLoader)
    # Flux.testmode!(model)
    all_predictions = []
    
    for data in testloader
        df, y,lat = data  # Assuming you don't have labels in the test data
        α, ŷ = model(df, :infer)
        save_predictions_to_nc(lat,α, ŷ,y, joinpath(@__DIR__,"predictions3_gc_oz.nc"))
        push!(all_predictions, ŷ)
    end
    
    return hcat(all_predictions...)
end

# Create an instance of the LinearHybridModel
model = LinearHybridModel_canopy(predictors, x, 1, 128)

# Train the model
checkpoint_file = "hybrid_clima_gc_oz.jld2"  # Specify the checkpoint file path here
train(model, trainloader, valloader, checkpoint_file);

load_model(model, checkpoint_file)
predictions = test_model(model, testloader) #valloader, testloader, trainall
