

cd("/Users/henriquemota/Library/CloudStorage/OneDrive-Personal/Documentos/_PhD_Classes/Trade Track/Spatial/Programming/")

using Pkg
Pkg.add(["CSV", "DataFrames"])

using CSV, DataFrames


n_counties = 3083



# Extracting the parameters for the counterfactual analysis 

#All these vectors are column vectors 
# Subscript n shall be row , Subscript i shall be second dim, col. 

D = CSV.read("output/D_n.csv", DataFrame)[:, 1]
R = CSV.read("output/R_n.csv", DataFrame)[:, 2] 
v = CSV.read("output/Vbar_n.csv", DataFrame)[:, 1]


w  =  CSV.read("output/w_i.csv", DataFrame)[:,1]
L = CSV.read("output/L_i.csv", DataFrame)[:, 2]
A = CSV.read("output/productivities.csv", DataFrame)[:, 1]

λ = CSV.read("output/lambda_ni.csv", DataFrame)[:,2:end]

hatB  =  ones(3083, 3083)
hat_κ =  ones(3083, 3083)

hat_v = function calc_Z(vn::Vector{Float64}, wage::Vector{Float64}, lambda::Matrix{Float64}, 
    hatB::Matrix{Float64}, hat_kappa::Vector{Float64},
    hat_w::Vector{Float64}, hat_lambda::Matrix{Float64})

sum_terms = hatB .* lambda .* ( reshape(wage, n_counties) ) .^ ϵ

hat_v = 

    return hatvn
    end 