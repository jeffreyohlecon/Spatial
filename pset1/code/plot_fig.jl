cd(project_path) 

default(dpi = 900)

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


# Calculate centiles for R
R = result_df_main.R
centiles_R = quantile(R, 0:0.01:1)
centiles_R[101] = centiles_R[101] + 1e-6

# Assign centiles to each county
result_df_main[!, :centile_R] = cut(R, centiles_R, labels=1:100)

# Calculate means of hat_L and hat_R for each centile of R
result_df_main[!, :adj_hatL] = (result_df_main[!, :hatL] .- 1) 
result_df_main[!, :adj_hatR] = (result_df_main[!, :hatR] .- 1) 

mean_hatL_by_centile_R = combine(groupby(result_df_main, :centile_R)) do df
    w = Weights(df.R)
    (; weighted_mean_adj_hatL = mean(df.adj_hatL, w))
end

mean_hatR_by_centile_R = combine(groupby(result_df_main, :centile_R)) do df
    w = Weights(df.R)
    (; weighted_mean_adj_hatR = mean(df.adj_hatR, w))
end

# Ensure centile_R is converted to numeric values
numeric_centile_R = parse.(Float64, string.(mean_hatR_by_centile_R.centile_R))


# Function to apply LOESS smoothing
function loess_smooth(x, y; span=0.5)
    model = loess(x, y; span=span)  # Fit LOESS model
    return predict(model, x)        # Get smoothed values
end

# Compute LOESS smoothed values
smooth_L = loess_smooth(numeric_centile_R, mean_hatL_by_centile_R.weighted_mean_adj_hatL)
smooth_R = loess_smooth(numeric_centile_R, mean_hatR_by_centile_R.weighted_mean_adj_hatR)


# Plot both curves in the same plot with different colors and include scatter points
plot(numeric_centile_R, smooth_L, label="", lw=2, color=:blue, xlabel="Centiles of R", ylabel="Weighted Mean Growth (%)", title="")
scatter!(numeric_centile_R, mean_hatL_by_centile_R.weighted_mean_adj_hatL, label="Labor", color=:blue, marker=:circle)
plot!(numeric_centile_R, smooth_R, label="", lw=2, color=:red)
scatter!(numeric_centile_R, mean_hatR_by_centile_R.weighted_mean_adj_hatR, label="Residents", color=:red, marker=:circle)

# Save the combined plot
savefig("output/figures/combined_scatter_hat_L_and_R_$version.png")

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
arlington_county_fips = 51013

# Get indices for the specified counties
cook_county_idx = findfirst(result_df_main.vector_number .== cook_county_fips)
nyc_county_idx = findfirst(result_df_main.vector_number .== nyc_county_fips)
san_diego_county_idx = findfirst(result_df_main.vector_number .== san_diego_county_fips)
san_bernardino_county_idx = findfirst(result_df_main.vector_number .== san_bernardino_fips)
middlesex_county_idx = findfirst(result_df_main.vector_number .== middlesex_county_fips) 
arlington_county_idx = findfirst(result_df_main.vector_number .== arlington_county_fips) 

# Add labels to the scatter plot


result_df_main.hatL[cook_county_idx]
result_df_main.hatR[cook_county_idx]

result_df_main

max_R_idx = argmax(result_df_main.R)

result_df_main.vector_number[max_R_idx]

xlims = (0.96, 1.08)
ylims = (0.96, 1.08)

# Fit a weighted linear model (Ordinary Least Squares)

weights = sqrt.(result_df_main.R)  # Assuming R contains the weights

X = [ones(length(result_df_main.hatL)) result_df_main.hatL]  # Design matrix
y = result_df_main.hatR  # Response variable

# Apply the weights to both X and y
X_weighted = X .* weights
y_weighted = y .* weights

# Solve the weighted least squares problem
weighted_coeffs = X_weighted \ y_weighted

# Compute fitted values
weighted_linear_fit = weighted_coeffs[1] .+ weighted_coeffs[2] .* result_df_main.hatL


scatter(result_df_main.hatL, result_df_main.hatR, markersize=sqrt.(R)/10^2, alpha=0.5,
    legend=false, xlabel="Hat L_i", ylabel="Hat R_i", 
    title="", xlims=xlims, ylims=ylims)
    # Add 45-degree dashed line
# Add horizontal and vertical lines at 1
hline!([1], color=:black, lw=2, linestyle=:dash, label="Horizontal line at 1")
vline!([1], color=:black, lw=2, linestyle=:dash, label="Vertical line at 1")
#Not sure I want this regression
#plot!(result_df_main.hatL, weighted_linear_fit, color=:blue, lw=2, label="Weighted Linear Fit")
# plot!([minimum(result_df_main.hatL), maximum(result_df_main.hatL)], 
    #    [minimum(result_df_main.hatL), maximum(result_df_main.hatL)], 
    #    color=:black, lw=2, linestyle=:dash, label="45-degree line")
        annotate!(result_df_main.hatL[cook_county_idx], result_df_main.hatR[cook_county_idx], text("Cook", :left, 8, :red))
        annotate!(result_df_main.hatL[nyc_county_idx], result_df_main.hatR[nyc_county_idx], text("New York", :left, 8, :red))
        annotate!(result_df_main.hatL[san_diego_county_idx], result_df_main.hatR[san_diego_county_idx], text("San Diego", :left, 8, :red))
      # annotate!(result_df_main.hatL[san_bernardino_county_idx], result_df_main.hatR[san_bernardino_county_idx], text("San Bernardino County", :left, 8, :red))
        annotate!(result_df_main.hatL[middlesex_county_idx], result_df_main.hatR[middlesex_county_idx], text("Middlesex", :left, 8, :red))
        annotate!(result_df_main.hatL[arlington_county_idx], result_df_main.hatR[arlington_county_idx], text("Arlington", :left, 8, :red))

# Save the plot
savefig("output/figures/scatterplot_hatL_hatR_categorized_bubble_sizes_$version.png")



function lorenz_curve(values)
    sorted_values = sort(values)
    cum_values = cumsum(sorted_values)
    total = sum(sorted_values)
    lorenz_points = cum_values ./ total
    return vcat(0.0, lorenz_points)  # Ensure correct shape
end




# Compute Lorenz curves
lorenz_hatR_R = lorenz_curve(Vector(result_df_main.hatR .* result_df_main.R .* result_df_main.hatv .* result_df_main.v))
lorenz_R = lorenz_curve(Vector(result_df_main.v .* result_df_main.R))

# Generate x-values
x_vals = range(0, stop=1, length=length(lorenz_hatR_R))

gini_after = 1 - 2 * sum(lorenz_hatR_R[2:end] .* diff(x_vals))
gini_before = 1 - 2 * sum(lorenz_R[2:end] .* diff(x_vals))

# Add text annotations for Gini indices
# Plot Lorenz curves
plot(x_vals, lorenz_hatR_R, label="Lorenz Curve of New Total Residential Income", lw=2, color=:blue, xlabel="Cumulative Share of Population", ylabel="Cumulative Share of Value", title="Lorenz Curve")
plot!(x_vals, lorenz_R, label="Lorenz Curve of Total Residential Income", lw=2, color=:red, linestyle=:dash)
plot!([0, 1], [0, 1], color=:black, lw=2, linestyle=:dash, label="45-degree line")
annotate!(0.15, 0.75, text("Gini Index Before: $(round(gini_before, digits=4))", :left, 10, :black))
annotate!(0.15, 0.65, text("Gini Index After: $(round(gini_after, digits=4))", :left, 10, :black))

# Save the plot
savefig("output/figures/lorenz_$version.png")

println("Generating other plots: $version")  
end

