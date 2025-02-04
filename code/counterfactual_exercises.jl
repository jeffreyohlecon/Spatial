using Pkg
using CSV, DataFrames, LinearAlgebra, Statistics, SparseArrays, Plots, StatsPlots, ColorSchemes, StatsBase
using CategoricalArrays

#Extracting some stuff from Dropbox before setting up the CD

#d = Matrix( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/d_ni.csv", DataFrame; header=false) )
hatB  =   Matrix( select( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/B_ni_hat_destination.csv", DataFrame), Not(:FIPS) ) )
π = Matrix(  select( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/pi_ni.csv", DataFrame), Not(:state_county_code) ) )
π = sparse(π)

cd("/Users/henriquemota/Library/CloudStorage/OneDrive-Personal/Documentos/_PhD_Classes/Trade Track/Spatial/Programming/")

param = CSV.read("output/commuting_parameters.csv", DataFrame)

ϵ = param[param.Parameter .== "epsilon", :Value][1]
ϕ = param[param.Parameter .== "phi", :Value][1]
α = param[param.Parameter .== "alpha", :Value][1]
σ = param[param.Parameter .== "sigma", :Value][1]
ψ = param[param.Parameter .== "psi", :Value][1]

# Extracting the parameters for the counterfactual analysis 

#All these vectors are column vectors 
# Subscript n shall be row , Subscript i shall be second dim, col. 

D =  Vector(CSV.read("output/D_n.csv", DataFrame, header = false)[:, 1])

L = Vector(CSV.read("output/L_i.csv", DataFrame, header = false)[:, 1])
R = Vector(CSV.read("output/R_n.csv", DataFrame, header = false)[:, 1])

w = Vector(CSV.read("output/w_i.csv", DataFrame, header = false)[:, 1])
v = Vector(CSV.read("output/Vbar_n.csv", DataFrame, header = false)[:, 1])

n_counties = size(L, 1)

Lbar = sum(L)

sum(L) - sum(R)

A = Vector(CSV.read("output/productivities.csv", DataFrame, header = false)[:, 1])

λ = Matrix( select( CSV.read("output/lambda_ni.csv", DataFrame), Not(:state_county_res) ) )

Y = w .* L 

# Shocks 
hat_κ =  ones(n_counties, n_counties)
hat_A =  ones(n_counties)
hat_d =  ones(n_counties, n_counties)

#Initial Guess for the Parameters 
guess_hat_wage = ones(n_counties)
guess_hat_lambda = ones(n_counties,n_counties)


function calc_hatv(vn::Vector{Float64}, wage::Vector{Float64}, lambda::Matrix{Float64}, #Observable 
    hatBni::Matrix{Float64}, hat_kappa::Matrix{Float64}, #Shock
    hat_wage::Vector{Float64}, hat_lambda::Matrix{Float64}) #Iterative Shock 

ek_term = hatBni .* lambda .* ( reshape(hat_wage, (1,n_counties)) ./ hat_kappa  ) .^ ϵ
# Mr. Eaton, Ora Pro Nobis 
denom_ek = sum(ek_term, dims = 2) 

sum_term = (ek_term ./ denom_ek ) .*  reshape(hat_wage, (1,n_counties)) .* reshape(wage, (1,n_counties)) 

x = (1 ./ vn) .* sum( sum_term, dims = 2 )

    return x[:,1]
    end 

function calc_LR(l::Vector{Int64}, r::Vector{Float64}, lambda::Matrix{Float64}, #Observable 
    hat_lambda::Matrix{Float64}) #Iterative Shock

    hatL  = ( sum(l) ./ reshape(l, (1, n_counties) ) ) .* sum( lambda .* hat_lambda, dims = 1)
    hatR = (sum(l) ./ r) .* sum( lambda .* hat_lambda, dims = 2)

    return hatL[1,:], hatR[:,1]
end 

function calc_hat_pi(p::SparseMatrixCSC{Float64},  #Observable
    hat_d::Matrix{Float64}, hat_A::Vector{Float64},  # Shock 
    hat_w::Vector{Float64}, hat_L::Vector{Float64}) #Iterative Shock
    
    ek_term = reshape(hat_L, (1, n_counties)) .* (hat_d .* reshape(hat_w./hat_A, (1, n_counties)) ) .^ (1-σ)
    denom = sum(p.*  ek_term, dims = 2)

   final = ek_term ./denom

return final
end

function calc_new_w(y::Vector{Float64}, p::SparseMatrixCSC{Float64}, v::Vector{Float64}, R::Vector{Float64},  D::Vector{Float64}, #Obsservable variables
   hat_pi::Matrix{Float64}, hat_L::Vector{Float64}, hat_R::Vector{Float64}, hat_v::Vector{Float64} ) #Iterative Shocks 

#This is how it looks in the files
num = sum(  p .* hat_pi .* ( hat_v .* hat_R .* v .* R .+ D), dims = 1) 
denom = reshape( y .* hat_L, (1, n_counties))    
new_w = num./denom

    return new_w 
end 


function calc_new_lambda(Y::Vector{Float64}, lambda::Matrix{Float64}, #Obsservable variables
    hatBni::Matrix{Float64}, hat_kappa::Matrix{Float64}, #Shock
    hat_P::Vector{Float64}, hat_Q::Vector{Float64}, hat_wage::Vector{Float64} ) #Iterative Shocks 
 
    QPterm =  ( (hat_P).^α  .* (hat_Q).^(1-α)  ) .^ (-ϵ) 
    w_kappa_term = (reshape(hat_wage, (1, n_counties) ) ./ hat_kappa) .^ (ϵ)

   num = hatBni .* QPterm .* w_kappa_term 
    
    denom = sum( lambda .* num)
 
    new_lambda = num ./ denom
 
    return new_lambda 
 end 

 ### NEW COUNTER-FACTUAL 

function loop_fun(old_wage::Vector{Float64}, old_lambda::Matrix{Float64}, #Past Wage, Lambda Guess (t+1) 
    hat_Bni::Matrix{Float64}, hat_kappa::Matrix{Float64}, hat_prod::Vector{Float64}, hat_dni::Matrix{Float64}, #Shocks
    bar_v::Vector{Float64},wage::Vector{Float64}, #Observable - Prices
    l::Vector{Int64}, r::Vector{Float64}, y::Vector{Float64}, deficit::Vector{Float64}, # Observables - Allocations
    lambda::Matrix{Float64}, p::SparseMatrixCSC{Float64, Int64} # Flow Data - Commuting and Trade 
)

hat_v = calc_hatv(bar_v, wage, lambda, hat_Bni, hat_kappa, old_wage, old_lambda)

hatL, hatR = calc_LR(l, r, lambda, old_lambda)

hatQ = hat_v .* hatR 
hat_pi = calc_hat_pi(p, hat_dni, hat_prod, old_wage, hatL)

hatP = ( reshape(hatL, (n_counties, 1)) ./ diag(hat_pi) ).^(1/(1-σ))   .* (old_wage ./ hat_prod)
hatP = hatP[:,1]

new_wage =  calc_new_w(y, p, v, r, deficit, hat_pi, hatL, hatR, hat_v) 
new_lambda = calc_new_lambda(y, lambda, hat_Bni, hat_kappa, hatP, hatQ, old_wage)

return reshape(new_wage, (n_counties, 1) ), new_lambda 
end 

# Shocks 
hat_κ =  ones(n_counties, n_counties)
hat_A =  ones(n_counties)
hat_d =  ones(n_counties, n_counties)

#Initial Guess for the Parameters 
guess_hat_wage = ones(n_counties)
guess_hat_lambda = ones(n_counties,n_counties)

#### Sanity Check: USE A SHOCK WITH ONES ONLY 
hat_B2 = ones(n_counties, n_counties) 

new_w, new_λ = loop_fun(guess_hat_wage, guess_hat_lambda, 
        hat_B2,  hat_κ, hat_A, hat_d, 
        v, w,
        L, R, Y, D, 
        λ, π)

new_w
new_λ

### START LOOP FOR COUNTERFACTUAL ANALYSIS 

tol = 1e-6
max_iter = 1000
iter = 0
ζ = 0.7

old_w = guess_hat_wage
old_λ = guess_hat_lambda

using Dates

start_time = now()

final_dist_w = 0.0
final_dist_λ = 0.0

while iter < max_iter

    tilde_w, tilde_λ = loop_fun(old_w, old_λ, 
        hatB, hat_κ, hat_A, hat_d, 
        v, w, L, R, Y, D, λ, π)
    
    dist_w = maximum( abs.(tilde_w ./ old_w .- 1 ) )
    dist_λ = maximum( abs.(tilde_λ ./ old_λ .- 1 ) )
    
    if dist_λ < tol && dist_w < tol
        final_dist_w = dist_w
        final_dist_λ = dist_λ
        break
    end
    
    old_w = tilde_w .* (1-ζ) + old_w .* ζ
    old_w = old_w[:,1]
    old_λ = tilde_λ .* (1-ζ) + old_λ .* ζ
    iter += 1
    if iter % 10 == 0
        println("Iteration $iter: max distance_w = $(maximum(abs.(dist_w))), max distance_λ = $(maximum(abs.(dist_λ)))")
    end
end

end_time = now()
elapsed_time = end_time - start_time

println("Converged after $iter iterations")
println("Final max distance_w = $final_dist_w, Final max distance_λ = $final_dist_λ")
println("Elapsed time: $elapsed_time")

old_w # I no these names don't make much sense
old_λ # I no these names don't make much sense

# Plot histogram and density for old_w
histogram(old_w, bins=100, normalize=true, alpha=0.5, label="Histogram of change in wages")
density!(old_w, label="Density of change in w")

# Plot histogram and density for old_λ
histogram(diag(old_λ), bins=100, normalize=true, alpha=0.5, label="Histogram of change in λ")
density!(old_λ[:], label="Density of λ")


cf_hatL, cf_hatR = calc_LR(L, R, λ, old_λ)

cf_hat_v = calc_hatv(v, w, λ, hatB, hat_κ, old_w, old_λ)

vector_number = select( CSV.read("output/lambda_ni.csv", DataFrame), :state_county_res)[:,1]

cf_hat_π = calc_hat_pi(π, hat_d, hat_A, old_w, cf_hatL)

cf_hatP = ( reshape(cf_hatL, (n_counties, 1)) ./ diag(cf_hat_π) ).^(1/(1-σ))   .* (old_w ./ hat_A)
cf_hatP = cf_hatP[:,1]

cf_hatQ = cf_hat_v .* cf_hatR 

real_v = cf_hat_v ./ cf_hatP



# Plot histogram and density for real_v
histogram(real_v, bins=100, normalize=true, alpha=0.5, label="Histogram of real wages")
density!(real_v, label="Density of real wages")

# Save the plot
savefig("output/figures/real_wages_histogram_density.png")


mean((old_w .* cf_hatL) ./ cf_hatP )
mean(old_w ./ cf_hatP)
mean(cf_hat_v ./ cf_hatP)
mean(cf_hatQ ./ cf_hatP)
mean(cf_hat_v ./ cf_hatP)

sum((cf_hatQ ./ cf_hatP .<1) ) # Half of counties worse off
# Calculate the weighted sum of real wages
sum((cf_hatQ ./ cf_hatP  .>1) .*  R ) ./sum(R) 
#Apparently very unequal, only 36% of original individuals are in counties that are better off


# Create a DataFrame with the required columns

# Preparing Dataframe for Plots
result_df = DataFrame(
    vector_number = vector_number,
    old_w = (old_w .-1) .*100, 
    cf_hatL = (cf_hatL .-1).*100, 
    cf_hatR = (cf_hatR .-1).*100, 
    diag_old_λ = (diag(old_λ) .-1).*100,
    real_v = (real_v .-1).*100
)

# Shocks for the maps 
CSV.write("output/map_database.csv", result_df)

result_df_main = DataFrame(
    vector_number = vector_number,
    hat_w = old_w,
    hatL = cf_hatL,
    hatR = cf_hatR,
    hatv = cf_hat_v, 
    diag_hat_λ = diag(old_λ), 
    hatP = cf_hatP,
    hatQ = cf_hatQ,
    real_v = real_v, 
    L = L, 
    R = R)

# Save the result_df as a .csv file in output as result_first_cf
CSV.write("output/result_first_cf.csv", result_df_main)

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

