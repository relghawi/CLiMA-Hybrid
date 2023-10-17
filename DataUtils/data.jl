using DataFrames
using AxisKeys
using DataFrameMacros
using Random
using Chain: @chain
using MLJ: partition
using Flux

function select_predictors(data, predictors)
    if isa(data, AbstractDataFrame)
        # If data is a DataFrame, select columns based on predictors
        return select(data, predictors) |> Matrix |> transpose
    elseif isa(data, DataFrameRow)
        # If data is a DataFrameRow, convert it to a DataFrame and then select columns
        df = DataFrame(data)
        return select(df, predictors) |> Matrix |> transpose
    elseif isa(data, AbstractVector)
        # If data is an AbstractVector, assume it's an array of data
        indices = [findfirst(isequal(predictor), predictors) for predictor in predictors]
        return data[indices]
    else
        throw(ArgumentError("Input data type not supported"))
    end
end

function select_predictors(df::KeyedArray, predictors)
    return df(predictors) |> Matrix
end

function select_cols(df::KeyedArray, predictors, x)
    # Merge predictors and x while ensuring uniqueness
    all_variables = unique([predictors..., x])

    # Select the desired columns from the KeyedArray
    selected_data = df(all_variables)

    return selected_data
end

# For a single variable
function select_variable(data, var::Symbol)
    if data isa KeyedArray
        return data(var) |> Vector
    elseif data isa DataFrame
        if var in propertynames(data)
            return data[!, var] |> Vector
        else
            throw(ArgumentError("Variable '$var' not found in the DataFrame"))
        end
    elseif isa(data, AbstractVector)
            return [data[findfirst(colname -> colname == var,vars)] for var in vars]
        end
end

# For multiple variables
function select_variable(data, vars::Vector{Symbol})
    if data isa KeyedArray
        return [data(var) |> Vector for var in vars]
    elseif data isa DataFrame
        return [data[!, var] |> Vector for var in vars]
    elseif isa(data, AbstractVector)
            return [data[findfirst(colname -> colname == var,vars)] for var in vars]
        else
            throw(ArgumentError("Variable(s) not found in the array of data"))
        end
end


function tokeyedArray(df)
    d = Matrix(df) |> transpose
    return KeyedArray(d, row=Symbol.(names(df)), col =1:size(d,2))
end

function select_cols(df, predictors, x)
    # Merge predictors and x while ensuring uniqueness
    all_variables = union(predictors, x)

    # Select the desired columns from the DataFrame
    selected_df = select(df, all_variables)

    return selected_df
end

## Create some data y = a ⋅ x1 + b, where a = f(x2,x3), b=2
# init dummy data
function gen_df(; seed = 123)
    Random.seed!(seed)
    df = DataFrame(rand(Float32, 1000, 3), :auto)
    # more variables
    df = @chain df begin
            @transform :a_syn = exp(-5f0((:x2-0.7f0))^2f0) + :x3/10f0
            @aside b = 2f0
            @transform :obs = :a_syn * :x1 + b
            @transform :pred_syn = :obs
            @transform :seqID = @bycol repeat(1:100, inner=10)
        end
     return df
end
df = gen_df()

function split_data(df, target, predictors, x; f = 0.8, batchsize=32, shuffle=true, partial=true)
    d_train, d_vali = partition(df, f; shuffle)
    # wrap training data into Flux.DataLoader
    df = select_cols(d_train, predictors, x)
    df = tokeyedArray(df)
    y = getproperty(d_train, target)
    data_t = (; df, y)
    trainloader = Flux.DataLoader(data_t; batchsize, shuffle, partial)
    trainall = Flux.DataLoader(data_t; batchsize=length(y), shuffle, partial)
    # wrap validation data into Flux.DataLoader
    df = select_cols(d_vali, predictors, x)
    df = tokeyedArray(df)
    y = getproperty(d_vali, target)
    data_v = (; df, y)
    valloader = Flux.DataLoader(data_v; batchsize=length(y), shuffle=false, partial=false)
    return trainloader, valloader, trainall, d_train, d_vali
end

function split_data(df, target, predictors, x; f = 0.8, v = 0.1, batchsize=32, shuffle=true, partial=true)
    d_train, d_temp = partition(df, f; shuffle)
    d_vali, d_test = partition(d_temp, v / (1 - f); shuffle)
    
    # Wrap training data into Flux.DataLoader
    df_train = select_cols(d_train, predictors, x)
    df_train = tokeyedArray(df_train)
    y_train = getproperty(d_train, target)
    data_train = (; df = df_train, y = y_train)
    trainloader = Flux.DataLoader(data_train; batchsize, shuffle, partial)
    trainall = Flux.DataLoader(data_train; batchsize=length(y_train), shuffle, partial)
    
    # Wrap validation data into Flux.DataLoader
    df_vali = select_cols(d_vali, predictors, x)
    df_vali = tokeyedArray(df_vali)
    y_vali = getproperty(d_vali, target)
    data_vali = (; df = df_vali, y = y_vali)
    valloader = Flux.DataLoader(data_vali; batchsize=length(y_vali), shuffle=false, partial=false)
    
    # Wrap test data into Flux.DataLoader
    df_test = select_cols(d_test, predictors, x)
    df_test = tokeyedArray(df_test)
    y_test = getproperty(d_test, target)
    data_test = (; df = df_test, y = y_test)
    testloader = Flux.DataLoader(data_test; batchsize=length(y_test), shuffle=false, partial=false)
    
    return trainloader, valloader, testloader, trainall, d_train, d_vali, d_test
end