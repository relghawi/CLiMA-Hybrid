module NormalizationModule
    export DataFrameNormalization, MatrixNormalization

    module DataFrameNormalization
        export DataFrameNormalizationStats, normalize, denormalize, fit_normalization, fit_and_normalize, fit_and_denormalize

        using DataFrames, Statistics
    
        struct DataFrameNormalizationStats
            mean::Vector{Float64}
            std::Vector{Float64}
        end
    
        function fit_normalization(data::DataFrame)
            mean_values = mean.(eachcol(data))
            std_values = std.(eachcol(data))
            stats = DataFrameNormalizationStats(mean_values, std_values)
            return stats
        end
    
        function normalize(stats::DataFrameNormalizationStats, data::DataFrame)
            num_cols = ncol(data)
            normalized_data = copy(data)  # Create a copy of the original DataFrame
        
            for col_idx in 1:num_cols
                # Get the mean value for the current column
                col_mean = stats.mean[col_idx]
        
                # Normalize each row in the current column using the corresponding mean
                normalized_data[:, col_idx] .= (data[:, col_idx] .- stats.mean[col_idx]) ./ stats.std[col_idx]
            end
                return normalized_data
        end
    
        function denormalize(stats::DataFrameNormalizationStats, data::DataFrame)
            num_cols = ncol(data)
            normalized_data = copy(data)  # Create a copy of the original DataFrame
        
            for col_idx in 1:num_cols
                # Get the mean value for the current column
                col_mean = stats.mean[col_idx]
        
                # Normalize each row in the current column using the corresponding mean
                denormalized_data[:, col_idx] .= (data[:, col_idx].* stats.std[col_idx]) .+ stats.mean[col_idx]
            end
                return denormalized_data
        end
    
    
        function fit_and_normalize(data::DataFrame)
            stats = fit_normalization(data)
            normalized_data = normalize(stats, data)
            return stats, normalized_data
        end
    
    
        function fit_and_denormalize(data::DataFrame)
        stats = fit_normalization(data)
        denormalized_data = denormalize(stats, data)
        return stats, denormalized_data
        end
    end

    module MatrixNormalizationModule
        export MatrixNormalizationStats, normalize, denormalize, fit_normalization, fit_and_normalize, fit_and_denormalize

        using DataFrames, Statistics

        struct MatrixNormalizationStats
            mean::Vector{Float64}
            std::Vector{Float64}
        end
    
        function fit_normalization(data::Matrix{Float64})
            mean_values = mean(data, dims=2)
            std_values = std(data, dims=2)
            stats = MatrixNormalizationStats(vec(mean_values), vec(std_values))
            return stats
        end
    
        function normalize(stats::MatrixNormalizationStats, data::Matrix{Float64})
            num_cols = size(data, 2)
            num_rows = size(data, 1)
            normalized_data = similar(data)
        
            for col_idx in 1:num_rows
                col_mean = stats.mean[col_idx]
                col_std = stats.std[col_idx]
        
                for row_idx in 1:num_cols
                    normalized_data[col_idx,row_idx] = (data[col_idx,row_idx] - col_mean) / col_std
                end
            end
    
            return normalized_data
        end
    
        function denormalize(stats::MatrixNormalizationStats,  data::Vector{Float32})
            num_cols = length(data)
            denormalized_data = similar(data)
                
            for col_idx in 1:num_cols

                denormalized_data[col_idx] = (data[col_idx] * stats.std[1]) + stats.mean[1]
            end
        
            return denormalized_data
        end
            
    
        function fit_and_normalize(data::Matrix{Float64})
            stats = fit_normalization(data)
            normalized_data = normalize(stats, data)
            return stats, normalized_data
        end
    
    end

end    
