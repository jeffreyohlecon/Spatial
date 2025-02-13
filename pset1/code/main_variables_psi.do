
	clear all
	set more off


	global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS1\CMLEE"
		global data "$path\data"
		global results "$path\results"
		global temp "$path\temp"

		
		insheet using "$data\data_for_reduced_form_fullTD.csv", clear
		sort state_county_code
		save "$temp\county_level_data.dta", replace

		
				
				
		********************************		
		* SECTION B.5 estimation of psi
		********************************

		use "$data\cfs2007.dta", clear

		*relation between the two measures of distance
		gen log_travelled_dist=log(travelled_dist)
		gen log_dist=log(dist)
		reg log_travelled_dist log_dist
		local sl=round(_b[log_dist],0.001)
		local r2=round(e(r2)*100,1)/100 //the simple rounding would add several zeros

		*filling in travelled distance based on relation with 
		*centroid distance for a few unreported pairs
		reg log_travelled_dist log_dist
		predict log_distance_hat if log_travelled_dist==.
		replace log_distance_hat=log_travelled_dist if log_travelled_dist~=.
		label var log_distance_hat "travelled distance + imputed when missing"


		*generating relevant variables
		gen log_value=log(value)
		**fixed effects
		egen orig=group(origin_cfs)
		egen dest=group(destination_cfs)

		*gravity slope is -1.29
		reghdfe log_value log_distance_hat, absorb(orig dest) vce(robust)
		
		* Extract estimated coefficient (β̂)
		local beta_hat = _b[log_distance_hat]
		
		* Given σ = 4, compute ψ
		local sigma = 4
		local psi = - `beta_hat' / (`sigma' - 1)

		* Save ψ in the CSV file
		file open results_file using "$results/estimated_parameters.csv", write append
		file write results_file "gravity_slope,`beta_hat'" _n
		file write results_file "psi,`psi'" _n
		file close results_file