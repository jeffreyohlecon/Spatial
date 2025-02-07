

# Calculate percentiles for L and R
percentiles_L = quantile(L, 0:0.05:1)
percentiles_R = quantile(R, 0:0.05:1)

percentiles_L[21] = percentiles_L[21] + 1e-6
percentiles_R[21] = percentiles_R[21] + 1e-6

# Assign percentiles to each county
result_df_main[!, :percentile_L] = cut(L, percentiles_L, labels=1:20)
result_df_main[!, :percentile_R] = cut(R, percentiles_R, labels=1:20)

# Calculate means of hat_L and hat_R for each percentile
result_df_main[!, :adj_hatL] = (result_df_main[!, :hatL] .- 1) .* 100
result_df_main[!, :adj_hatR] = (result_df_main[!, :hatR] .- 1) .* 100

mean_hatL_by_percentile_L = combine(groupby(result_df_main, :percentile_L)) do df
    w = Weights(df.L)
    (; weighted_mean_adj_hatL = mean(df.adj_hatL, w))
end

# Compute weighted mean of adj_hatR using R as weight
mean_hatR_by_percentile_R = combine(groupby(result_df_main, :percentile_R)) do df
    w = Weights(df.R)
    (; weighted_mean_adj_hatR = mean(df.adj_hatR, w))
end

# Generate bar plots for weighted mean_hatL based on percentiles of L
bar_plot_weighted_L = @df mean_hatL_by_percentile_L bar(:percentile_L, :weighted_mean_adj_hatL, legend=false, xlabel="Percentiles of L", ylabel="Weighted Mean Growth", title="Weighted Mean Growth by Percentiles of L")

# Generate bar plots for weighted mean_hatR based on percentiles of R
bar_plot_weighted_R = @df mean_hatR_by_percentile_R bar(:percentile_R, :weighted_mean_adj_hatR, legend=false, xlabel="Percentiles of R", ylabel="Weighted Mean Growth", title="Weighted Mean Growth by Percentiles of R")

# Save the plots
savefig(bar_plot_weighted_L, "output/figures/weighted_hat_L_by_percentiles_of_L.png")
savefig(bar_plot_weighted_R, "output/figures/weighted_hat_R_by_percentiles_of_R.png")

# Calculate deciles for L and R
deciles_L = quantile(L, 0:0.1:1)
deciles_R = quantile(R, 0:0.1:1)

deciles_L[11] = deciles_L[11] + 1e-6
deciles_R[11] = deciles_R[11] + 1e-6

# Assign deciles to each county
result_df_main[!, :decile_L] = cut(L, deciles_L, labels=1:10)
result_df_main[!, :decile_R] = cut(R, deciles_R, labels=1:10)

# Calculate means of hat_L and hat_R for each decile
result_df_main[!, :adj_hatL] = (result_df_main[!, :hatL] .- 1) .* 100
result_df_main[!, :adj_hatR] = (result_df_main[!, :hatR] .- 1) .* 100

# Calculate means of hat_L and hat_R for each decile
mean_hatL_by_decile_L = combine(groupby(result_df_main, :decile_L)) do df
    w = Weights(df.L)
    (; weighted_mean_adj_hatL = mean(df.adj_hatL, w))
end

mean_hatR_by_decile_R = combine(groupby(result_df_main, :decile_R)) do df
    w = Weights(df.R)
    (; weighted_mean_adj_hatR = mean(df.adj_hatR, w))
end

# Generate bar plots for weighted mean_hatL based on deciles of L
bar_plot_weighted_L_decile = @df mean_hatL_by_decile_L bar(:decile_L, :weighted_mean_adj_hatL, legend=false, xlabel="Deciles of L", ylabel="Weighted Mean Growth", title="Weighted Mean Growth by Deciles of L")

# Generate bar plots for weighted mean_hatR based on deciles of R
bar_plot_weighted_R_decile = @df mean_hatR_by_decile_R bar(:decile_R, :weighted_mean_adj_hatR, legend=false, xlabel="Deciles of R", ylabel="Weighted Mean Growth", title="Weighted Mean Growth by Deciles of R")

# Save the plots
savefig(bar_plot_weighted_L_decile, "output/figures/weighted_hat_L_by_deciles_of_L.png")
savefig(bar_plot_weighted_R_decile, "output/figures/weighted_hat_R_by_deciles_of_R.png")