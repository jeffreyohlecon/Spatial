using Pkg
using CSV, DataFrames, LinearAlgebra, Statistics, SparseArrays, Plots, StatsPlots, ColorSchemes, StatsBase, GLM
using Dates
using CategoricalArrays
using Loess 
using DataFrames, StatsPlots, StatsBase, CategoricalArrays
using Shapefile
using DataFrames
using Plots
using CSV 
import .Missings 

project_path = joinpath(homedir(), "Library", "CloudStorage", "OneDrive-Personal", 
                        "Documentos", "_PhD_Classes", "Trade Track", "Spatial", "Programming")

# Stop execution if the input is empty
isempty(project_path) && error("You must enter a valid directory path!")

# Define the versions to loop over
versions = ["uniform", "het"]
ones_diagonal = "yes" # "yes" or "no"

for version in versions    
println("Running counterfactual analysis for version: $version")

#Extracting some stuff from Dropbox before setting up the CD

#d = Matrix( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/d_ni.csv", DataFrame; header=false) )
dropbox_path = joinpath(homedir(), "Dropbox", "BigDataFiles") 

cd(dropbox_path)

if version == "uniform"
    hatB  =   Matrix( select( CSV.read("B_ni_hat_uniform.csv", DataFrame), Not(:FIPS) ) )
end
if version == "het"
hatB  =   Matrix( select( CSV.read("B_ni_hat_destination_sorted.csv", DataFrame), Not(:FIPS) ) )
end

if ones_diagonal == "yes"
println("ONES ARE IN THE DIAGONAL")
for i in 1:size(hatB, 1)
    hatB[i, i] = 1
end
else 
println("ONES ARE NOT IN THE DIAGONAL")
end 

π = Matrix(  select( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/pi_ni.csv", DataFrame), Not(:state_county_code) ) )
π = sparse(π)

cd(project_path) 

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

san_check1 = new_w[1:10] 
san_check2 = new_λ[1:5, 1:5]
println("Sanity Check - Shock of 1: wage (First 10 entries): $san_check1")
println("Sanity Check - Shock of 1: lambda (First 5x5 entries): $san_check2")


### START LOOP FOR COUNTERFACTUAL ANALYSIS 

tol = 1e-6
max_iter = 1000
ζ = 0.7

function eq_gen(guess_w::Vector{Float64}, guess_lambda::Matrix{Float64}, tol::Float64, max_iter::Int64, ζ::Float64,
    hatB::Matrix{Float64}, hat_κ::Matrix{Float64}, hat_A::Vector{Float64}, hat_d::Matrix{Float64}, 
    v::Vector{Float64}, w::Vector{Float64}, L::Vector{Int64}, R::Vector{Float64}, Y::Vector{Float64}, D::Vector{Float64},
     λ::Matrix{Float64}, π::SparseMatrixCSC{Float64, Int64})

start_time = now()

iter = 0
final_dist_w = 0 
final_dist_λ = 0
old_w = guess_w 
old_λ = guess_lambda

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
return old_w, old_λ
end 

cf_w, cf_λ = eq_gen(guess_hat_wage, guess_hat_lambda, tol, max_iter, ζ, hatB, hat_κ, hat_A, hat_d, v, w, L, R, Y, D, λ, π)
cf_hatL, cf_hatR = calc_LR(L, R, λ, cf_λ)

cf_hat_v = calc_hatv(v, w, λ, hatB, hat_κ, cf_w, cf_λ)

vector_number = select( CSV.read("output/lambda_ni.csv", DataFrame), :state_county_res)[:,1]

cf_hat_π = calc_hat_pi(π, hat_d, hat_A, cf_w, cf_hatL)

cf_hatP = ( reshape(cf_hatL, (n_counties, 1)) ./ diag(cf_hat_π) ).^(1/(1-σ))   .* (cf_w ./ hat_A)
cf_hatP = cf_hatP[:,1]

cf_hatQ = cf_hat_v .* cf_hatR 

real_v = cf_hat_v ./ cf_hatP

println("Sanity Check Results:")
println("Sum of cf_λ * λ: ", sum(cf_λ .* λ))
println("Sum of cf_hat_π * π: ", sum(cf_hat_π .* π))
println("Sum of cf_hatL * L / sum(L): ", sum(cf_hatL .* L) / sum(L))
println("Sum of cf_hatR * R / sum(R): ", sum(cf_hatR .* R) / sum(R))
#
var_hatB = sum( (hatB).^(1/ϵ) .* λ, dims =2) ./ sum(λ, dims = 2)

non_demean_shock = sum( (hatB).^(1/ϵ) .* λ, dims =2) ./ sum(λ, dims = 2)

non_demean_shock = non_demean_shock[:,1]

var_hatB = ( var_hatB .- mean(var_hatB))[:,1]

# Create a DataFrame with the required columns

# Preparing Dataframe for Plots
result_df = DataFrame(
    vector_number = vector_number,
    cf_w = (cf_w .-1) .*100, 
    cf_hatL = (cf_hatL .-1).*100, 
    cf_hatR = (cf_hatR .-1).*100, 
    diag_cf_λ = (diag(cf_λ) .-1).*100,
    real_v = (real_v .-1).*100,
    shock = var_hatB,
    non_demean_shock = non_demean_shock 
)

# Shocks for the maps 
CSV.write("output/map_database_$version.csv", result_df)

result_df_main = DataFrame(
    vector_number = vector_number,
    hat_w = cf_w,
    hatL = cf_hatL,
    hatR = cf_hatR,
    hatv = cf_hat_v, 
    diag_hat_λ = diag(cf_λ), 
    diag_hat_pi = diag(cf_hat_π), 
    hatP = cf_hatP,
    hatQ = cf_hatQ,
    real_v = real_v, 
    var_hatB = var_hatB ,   
    L = L, 
    R = R,
    A = A, 
    w = w, 
    non_demean_shock = non_demean_shock  
    )

# Save the result_df as a .csv file in output as result_first_cf
CSV.write("output/result_cf_$version.csv", result_df_main)

U = ( (1 ./diag(cf_λ)).^(1/ϵ) ) .* ( (1 ./ diag(cf_hat_π)).^(α/(σ-1)) ).* 
((cf_w ./ cf_hat_v).^ (1-α)).*((cf_hatL).^(α/(σ-1)) .* (cf_hatR).^(α-1) )

welfare_gain  = maximum(U)
println("Welfare Gain: $welfare_gain")

# Save the welfare gain to a CSV file
welfare_gain_df = DataFrame(welfare_gain = welfare_gain)
CSV.write("output/welfare_gain_$version.csv", welfare_gain_df)

# Plot histogram and density for U
histogram(U, bins=100, normalize=true, alpha=0.5, label="Histogram of U")
density!(U, label="Density of U")

# Save the plot
savefig("output/figures/histogram_density_U_$version.png") 

println("Counter-factual analysis finished for version: $version")  
end

# include("plot_maps.jl")
# include("plot_fig.jl")

