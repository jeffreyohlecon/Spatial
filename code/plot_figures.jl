# Pkg.add(["Shapefile", "DataFrames", "FilePathsBase", "Plots", "CSV", "Missings"])

using Shapefile
using DataFrames
using Plots
using CSV 
import .Missings 

# Define the shapefile path in a way that is adaptable to any computer
shapefile_path = joinpath(homedir(), "Dropbox", "BigDataFiles", "2010_county_shape_file", "tl_2010_us_county10.shp") 

# Define your Directory here 
cd("/Users/henriquemota/Library/CloudStorage/OneDrive-Personal/Documentos/_PhD_Classes/Trade Track/Spatial/Programming/")

# Load the CSV file into a DataFrame
result_df = DataFrame(CSV.File("output/map_database.csv"))
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

    p = plot(color=:blues, grid=false, framestyle=:none, axis=nothing, ratio=1, 
             xlims=(xmin, xmax), ylims=(ymin, ymax), title=title_name, size = (900,675))

            
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
p1 = plot_variable(merged_df, :old_w, "Old Wages by County")
p2 = plot_variable(merged_df, :cf_hatL, "Change in Workers by County")
p3 = plot_variable(merged_df, :cf_hatR, "Change in Residents by County")
p4 = plot_variable(merged_df, :diag_old_λ, "Change in Living Where You Work by County")
p5 = plot_variable(merged_df, :real_v, "Change in Real Average Income by County")

# Save the plots to the output/figures folder with higher resolution
savefig(p1, "output/figures/old_wages_by_county.png")
savefig(p2, "output/figures/change_in_workers_by_county.png")
savefig(p3, "output/figures/change_in_residents_by_county.png")
savefig(p4, "output/figures/change_in_living_where_you_work_by_county.png")
savefig(p5, "output/figures/change_in_real_income.png")

