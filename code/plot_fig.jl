cd(project_path) 

versions = ["uniform", "het"]

for version in versions    
println("Running figures for version: $version")

result_df_main =  DataFrame( CSV.File("output/result_cf_$version.csv") )

real_v = result_df_main.real_v

# Plot histogram and density for real_v
histogram(result_df_main.real_v, bins=100, normalize=true, alpha=0.5, label="Histogram of Real Wages")
density!(real_v, label="Density of real wages")

# Save the plot
savefig("output/figures/real_wages_histogram_density_$version.png")

# Calculate centiles for L and R
L = result_df_main.L
R = result_df_main.R

centiles_L = quantile(L, 0:0.01:1)
centiles_R = quantile(R, 0:0.01:1)

centiles_L[101] = centiles_L[101] + 1e-6
centiles_R[101] = centiles_R[101] + 1e-6

# Assign centiles to each county
result_df_main[!, :centile_L] = cut(L, centiles_L, labels=1:100)
result_df_main[!, :centile_R] = cut(R, centiles_R, labels=1:100)

# Calculate means of hat_L and hat_R for each centile
result_df_main[!, :adj_hatL] = (result_df_main[!, :hatL] .- 1) .* 100
result_df_main[!, :adj_hatR] = (result_df_main[!, :hatR] .- 1) .* 100

mean_hatL_by_centile_L = combine(groupby(result_df_main, :centile_L)) do df
    w = Weights(df.L)
    (; weighted_mean_adj_hatL = mean(df.adj_hatL, w))
end

mean_hatR_by_centile_R = combine(groupby(result_df_main, :centile_R)) do df
    w = Weights(df.R)
    (; weighted_mean_adj_hatR = mean(df.adj_hatR, w))
    end

# Ensure centile_L is converted to numeric values
numeric_centile_L = parse.(Float64, string.(mean_hatL_by_centile_L.centile_L))
numeric_centile_R = parse.(Float64, string.(mean_hatR_by_centile_R.centile_R))

# Function to apply LOESS smoothing
function loess_smooth(x, y; span=0.5)
    model = loess(x, y; span=span)  # Fit LOESS model
    return predict(model, x)        # Get smoothed values
end

# Compute LOESS smoothed values
smooth_L = loess_smooth(numeric_centile_L, mean_hatL_by_centile_L.weighted_mean_adj_hatL)
smooth_R = loess_smooth(numeric_centile_R, mean_hatR_by_centile_R.weighted_mean_adj_hatR)

# Scatter plot with LOESS smoothing for L
scatter_plot_weighted_L = scatter(numeric_centile_L, mean_hatL_by_centile_L.weighted_mean_adj_hatL, 
    legend=false, xlabel="Centiles of L", ylabel="Weighted Mean Growth (%)", 
    title="Weighted Mean Growth (%) by Centiles of L", 
    xticks=([0, 10, 25, 50, 75, 90, 100], string.([0, 10, 25, 50, 75, 90,100])))

plot!(numeric_centile_L, smooth_L, label="loess", lw=2, color=:blue)  # Add LOESS line

# Scatter plot with LOESS smoothing for R
scatter_plot_weighted_R = scatter(numeric_centile_R, mean_hatR_by_centile_R.weighted_mean_adj_hatR, 
    legend=false, xlabel="Centiles of R", ylabel="Weighted Mean Growth (%)", 
    title="Weighted Mean Growth (%) by Centiles of R", 
    xticks=([0, 10, 25, 50, 75, 90, 100], string.([0, 10, 25, 50, 75, 90,100])))

plot!(numeric_centile_R, smooth_R, label="loess", lw=2, color=:blue)

    # Save the plots
savefig(scatter_plot_weighted_L, "output/figures/scatter_hat_L_$version.png")
savefig(scatter_plot_weighted_R, "output/figures/scatter_hat_R_$version.png")

# Filter out high outliers in A
A = result_df_main.A
var_hatB =  result_df_main.var_hatB

filtered_indices = abs.(A .- mean(A)) .< 3 * std(A)
filtered_var_hatB = var_hatB[filtered_indices]
filtered_A = A[filtered_indices]

# Fit a linear model (Ordinary Least Squares)
coeffs = [ones(length(filtered_A)) filtered_A] \ filtered_var_hatB  # Solve Ax = b
linear_fit = coeffs[1] .+ coeffs[2] .* filtered_A  # Compute fitted values

# Scatter plot with linear fit
scatter(filtered_A, filtered_var_hatB, legend=false, xlabel="A", ylabel="var_hatB", 
    title="Scatter plot of var_hatB on A with Linear Fit")

plot!(filtered_A, linear_fit, color=:red, lw=2, label="")  # Add linear fit line

savefig("output/figures/scatterplot_with_linear_fit_$version.png")

# Calculate deciles for R
#pctile_R = quantile(R, 0:0.25:1)
#pctile_R[5] = pctile_R[5] + 1e-6

# Assign deciles to each county
#result_df_main[!, :pctile_R] = cut(R, pctile_R, labels=1:4)

# Update centile_R_values to decile_R values
#pctile_val = convert(Vector{Float64}, result_df_main.pctile_R)

# Generate scatter plot of hat_L and hat_R with categorized bubble sizes of R


# Define FIPS codes for the counties
cook_county_fips = 17031
nyc_county_fips = 36061
san_diego_county_fips = 6073
san_bernardino_fips = 6071
middlesex_county_fips = 25017

# Get indices for the specified counties
cook_county_idx = findfirst(result_df_main.vector_number .== cook_county_fips)
nyc_county_idx = findfirst(result_df_main.vector_number .== nyc_county_fips)
san_diego_county_idx = findfirst(result_df_main.vector_number .== san_diego_county_fips)
san_bernardino_county_idx = findfirst(result_df_main.vector_number .== san_bernardino_fips)
middlesex_county_idx = findfirst(result_df_main.vector_number .== middlesex_county_fips) 

# Add labels to the scatter plot


result_df_main.hatL[cook_county_idx]
result_df_main.hatR[cook_county_idx]

result_df_main

max_R_idx = argmax(result_df_main.R)

result_df_main.vector_number[max_R_idx]

scatter(result_df_main.hatL, result_df_main.hatR, markersize=sqrt.(R)/10^2, alpha=0.5,
    legend=false, xlabel=" Hat L_i", ylabel="Hat R_i", 
    title="")
    # Add 45-degree dashed line

    plot!([minimum(result_df_main.hatL), maximum(result_df_main.hatL)], 
        [minimum(result_df_main.hatL), maximum(result_df_main.hatL)], 
        color=:black, lw=2, linestyle=:dash, label="45-degree line")
        annotate!(result_df_main.hatL[cook_county_idx], result_df_main.hatR[cook_county_idx], text("Cook County", :left, 8, :red))
        annotate!(result_df_main.hatL[nyc_county_idx], result_df_main.hatR[nyc_county_idx], text("NYC County", :left, 8, :red))
        annotate!(result_df_main.hatL[san_diego_county_idx], result_df_main.hatR[san_diego_county_idx], text("San Diego County", :left, 8, :red))
      # annotate!(result_df_main.hatL[san_bernardino_county_idx], result_df_main.hatR[san_bernardino_county_idx], text("San Bernardino County", :left, 8, :red))
        annotate!(result_df_main.hatL[middlesex_county_idx], result_df_main.hatR[middlesex_county_idx], text("Middlesex County", :left, 8, :red))


# Save the plot
savefig("output/figures/scatterplot_hatL_hatR_categorized_bubble_sizes_$version.png")

println("Generating other plots: $version")  
end
