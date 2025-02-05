# Define the shapefile path in a way that is adaptable to any computer
shapefile_path = joinpath(homedir(), "Dropbox", "BigDataFiles", "2010_county_shape_file", "tl_2010_us_county10.shp") 

# Define your Directory here 
cd(project_path)

versions = ["uniform", "het"]

for version in versions    
    println("Running figures for version: $version")

# Load the CSV file into a DataFrame
result_df = DataFrame(CSV.File("output/map_database_$version.csv"))
shapes = Shapefile.Table(shapefile_path)

# Convert the shapefile to a DataFrame
shape_df = DataFrame(shapes)

# Ensure vector_number has 5 digits by adding leading zeros if necessary
# Create a new vector with 5-digit FIPS codes
result_df.vector_number_5digit = lpad.(string.(result_df.vector_number), 5, '0')

# Merge the result_df with the shape_df on the FIPS code
merged_df = leftjoin(shape_df, result_df, on = :GEOID10 => :vector_number_5digit)

default(dpi = 900)

# Define geographic bounds for contiguous U.S.
xmin, xmax = -125, -65  
ymin, ymax = 24, 50      

# Function to generate each map (returns a plot object)
function plot_variable(data, variable_name, title_name)

    p = plot(color=:blues, grid=false, framestyle=:box, axis=nothing, ratio=1, 
             xlims=(xmin, xmax), ylims=(ymin, ymax), title=title_name)

            
min_val = minimum(skipmissing(data[!,variable_name]))
max_val = maximum(skipmissing(data[!,variable_name]))

val = maximum([abs(min_val), abs(max_val)]) 

clims=(-1.05, 1.05).*val

# Define the custom 3-color scale with min_val, 0, and max_val
custom_cmap = cgrad([:red, :white, :green], scale=:symlog, rev=false)

    for row in eachrow(data)
        if !ismissing(row[variable_name])
           # plot!(p, row.geometry, fill_z=row[variable_name], label="", lw=0.5, linecolor=:black, c=custom_cmap, clims=clims, colorbar_size=0.03, colorbar_title="", colorbar_titlefont=font(8))
            plot!(p, row.geometry, fill_z=row[variable_name], label="", lw=0.5, linecolor=:black, c=custom_cmap, clims = clims, colorbar_size=0.03)
        end
    end

    return p  # Ensure function returns a plot object
end

# Generate maps for the four variables
# Use the PlotlyJS backend for interactive plots
p1 = plot_variable(merged_df, :cf_w, "Change in Wages (%) by County")
p2 = plot_variable(merged_df, :cf_hatL, "Change in Workers (%) by County")
p3 = plot_variable(merged_df, :cf_hatR, "Change in Residents (%) by County")
p4 = plot_variable(merged_df, :diag_cf_λ, "Change in Living Where You Work (%) by County")
p5 = plot_variable(merged_df, :real_v, "Change in Real Average Income (%) by County")
p6 = plot_variable(merged_df, :shock, "Change in Average Welfare - Hat Bar Un")

# Save the plots to the output/figures folder with higher resolution
savefig(p1, "output/figures/cf_wages_by_county_$version.png")
savefig(p2, "output/figures/change_in_workers_by_county_$version.png")
savefig(p3, "output/figures/change_in_residents_by_county_$version.png")
savefig(p4, "output/figures/change_in_living_where_you_work_by_county_$version.png")
savefig(p5, "output/figures/change_in_real_income_$version.png")
savefig(p6, "output/figures/change_in_B_$version.png")


# Filter the DataFrame for New York state (FIPS code for New York is '36')
ny_df = filter(row -> startswith(row.GEOID10, "36"), merged_df)

# Define geographic bounds for New York
xmin_ny, xmax_ny = -80, -70  
ymin_ny, ymax_ny = 40, 45      

# Function to generate each map for New York (returns a plot object)
function plot_variable_ny(data, variable_name, title_name)

    p = plot(color=:blues, grid=false, framestyle=:box, axis=nothing, ratio=1, 
             xlims=(xmin_ny, xmax_ny), ylims=(ymin_ny, ymax_ny), title=title_name)

    min_val = minimum(skipmissing(data[!,variable_name]))
    max_val = maximum(skipmissing(data[!,variable_name]))

    val = maximum([abs(min_val), abs(max_val)]) 

    clims=(-1.05, 1.05).*val

    # Define the custom 3-color scale with min_val, 0, and max_val
    custom_cmap = cgrad([:red, :white, :green], scale=:symlog, rev=false)

    for row in eachrow(data)
        if !ismissing(row[variable_name])
            plot!(p, row.geometry, fill_z=row[variable_name], label="", lw=0.5, linecolor=:black, c=custom_cmap, clims = clims, colorbar_size=0.03)
        end
    end

    return p  # Ensure function returns a plot object
end

# Generate maps for the four variables for New York
p1_ny = plot_variable_ny(ny_df, :cf_w, "Change in Wages (%) by County")
p2_ny = plot_variable_ny(ny_df, :cf_hatL, "Change in Workers (%) by County")
p3_ny = plot_variable_ny(ny_df, :cf_hatR, "Change in Residents (%) by County")
p4_ny = plot_variable_ny(ny_df, :diag_cf_λ, "Change in Living Where You Work (%) by County")
p5_ny = plot_variable_ny(ny_df, :real_v, "Change in Real Average Income (%) by County")
p6_ny = plot_variable_ny(ny_df, :shock, "Change in Average Welfare - Hat Bar Un")

# Save the plots to the output/figures folder with higher resolution
savefig(p1_ny, "output/figures/new_york/ny_cf_wages_by_county_$version.png")
savefig(p2_ny, "output/figures/new_york/ny_change_in_workers_by_county_$version.png")
savefig(p3_ny, "output/figures/new_york/ny_change_in_residents_by_county_$version.png")
savefig(p4_ny, "output/figures/new_york/ny_change_in_living_where_you_work_by_county_$version.png")
savefig(p5_ny, "output/figures/new_york/ny_change_in_real_income_$version.png")
savefig(p6_ny, "output/figures/new_york/ny_change_in_shock_$version.png")

# Filter the DataFrame for California state (FIPS code for California is '06')
ca_df = filter(row -> startswith(row.GEOID10, "06"), merged_df)

# Define geographic bounds for California
xmin_ca, xmax_ca = -125, -114  
ymin_ca, ymax_ca = 32, 42      

# Function to generate each map for California (returns a plot object)
function plot_variable_ca(data, variable_name, title_name)

    p = plot(color=:blues, grid=false, framestyle=:box, axis=nothing, ratio=1, 
             xlims=(xmin_ca, xmax_ca), ylims=(ymin_ca, ymax_ca), title=title_name)

    min_val = minimum(skipmissing(data[!,variable_name]))
    max_val = maximum(skipmissing(data[!,variable_name]))

    val = maximum([abs(min_val), abs(max_val)]) 

    clims=(-1.05, 1.05).*val

    # Define the custom 3-color scale with min_val, 0, and max_val
    custom_cmap = cgrad([:red, :white, :green], scale=:symlog, rev=false)

    for row in eachrow(data)
        if !ismissing(row[variable_name])
            plot!(p, row.geometry, fill_z=row[variable_name], label="", lw=0.5, linecolor=:black, c=custom_cmap, clims = clims, colorbar_size=0.03)
        end
    end

    return p  # Ensure function returns a plot object
end

# Generate maps for the four variables for California
p1_ca = plot_variable_ca(ca_df, :cf_w, "Change in Wages (%) by County")
p2_ca = plot_variable_ca(ca_df, :cf_hatL, "Change in Workers (%) by County")
p3_ca = plot_variable_ca(ca_df, :cf_hatR, "Change in Residents (%) by County")
p4_ca = plot_variable_ca(ca_df, :diag_cf_λ, "Change in Living Where You Work (%) by County")
p5_ca = plot_variable_ca(ca_df, :real_v, "Change in Real Average Income (%) by County")
p6_ca = plot_variable_ca(ca_df, :shock, "Change in Average Welfare - Hat Bar Un")

# Save the plots to the output/figures folder with higher resolution
savefig(p1_ca, "output/figures/california/ca_cf_wages_by_county_$version.png")
savefig(p2_ca, "output/figures/california/ca_change_in_workers_by_county_$version.png")
savefig(p3_ca, "output/figures/california/ca_change_in_residents_by_county_$version.png")
savefig(p4_ca, "output/figures/california/ca_change_in_living_where_you_work_by_county_$version.png")
savefig(p5_ca, "output/figures/california/ca_change_in_real_income_$version.png")
savefig(p6_ca, "output/figures/california/ca_change_in_B_$version.png")

end
