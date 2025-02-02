using Pkg
Pkg.add(["CSV", "DataFrames"])
using CSV, DataFrames, LinearAlgebra, Statistics, SparseArrays

#Extracting some stuff from Dropbox before setting up the CD

#d = Matrix( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/d_ni.csv", DataFrame; header=false) )
hatB  =   Matrix( select( CSV.read("/Users/henriquemota/Dropbox/BigDataFiles/B_ni_hat.csv", DataFrame), Not(:FIPS) ) )
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

D = CSV.read("output/D_n.csv", DataFrame)[!, :D_n].* 10^6
# D = CSV.read("output/data_esteban_rf.csv", DataFrame)[!, :deficit]


esteban = CSV.read("output/data_esteban_rf.csv", DataFrame)


R = Vector(CSV.read("output/R_n.csv", DataFrame, header = false)[:, 1])
v = Vector(CSV.read("output/Vbar_n.csv", DataFrame, header = false)[:, 1])

w = Vector(CSV.read("output/w_i.csv", DataFrame, header = false)[:, 1])
L = Vector(CSV.read("output/L_i.csv", DataFrame, header = false)[:, 1])

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




function calc_new_w(y::Vector{Float64}, p::SparseMatrixCSC{Float64}, v::Vector{Float64}, R::Vector{Float64}, #Obsservable variables
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





# Shocks 

hat_κ =  ones(n_counties, n_counties)
hat_A =  ones(n_counties)
hat_d =  ones(n_counties, n_counties)

#Initial Guess for the Parameters 
guess_hat_wage = ones(n_counties)
guess_hat_lambda = ones(n_counties,n_counties)

#### Sanity Check: USE A SHOCK WITH ONES ONLY 
hat_B2 = ones(n_counties, n_counties)

hat_v = calc_hatv(v, w, λ, hat_B2, hat_κ, guess_hat_wage, guess_hat_lambda)

hatL, hatR = calc_LR(L, R, λ, guess_hat_lambda)

hatQ = hat_v .* hatR 

hat_pi = calc_hat_pi(π, hat_d, hat_A, guess_hat_wage, hatL)

hatP = ( reshape(hatL, (n_counties, 1)) ./ diag(hat_pi) ).^(1/(1-σ))   .* (guess_hat_wage ./ hat_A)
hatP = hatP[:,1]

new_wage =  calc_new_w(Y, π, v, R, hat_pi, hatL, hatR, hat_v) 
new_lambda = calc_new_lambda(Y, λ, hat_B2, hat_κ, hatP, hatQ, guess_hat_wage)





