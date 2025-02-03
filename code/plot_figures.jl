using Shapefile
using DataFrames
using Plots

# Preparing Dataframe for Plots
result_df = DataFrame(
    vector_number = vector_number,
    old_w = (old_w .-1) .*100, 
    cf_hatL = (cf_hatL .-1).*100, 
    cf_hatR = (cf_hatR .-1).*100, 
    diag_old_λ = (diag(old_λ) .-1).*100
)

# Load the shapefile
shapefile_path = "/Users/henriquemota/Dropbox/BigDataFiles/2010_county_shape_file/tl_2010_us_county10.shp"
shapes = Shapefile.Table(shapefile_path)

# Convert the shapefile to a DataFrame
shape_df = DataFrame(shapes)

# Ensure vector_number has 5 digits by adding leading zeros if necessary
# Create a new vector with 5-digit FIPS codes
result_df.vector_number_5digit = lpad.(string.(vector_number), 5, '0')

# Merge the result_df with the shape_df on the FIPS code
merged_df = leftjoin(shape_df, result_df, on = :GEOID10 => :vector_number_5digit)

using Plots

default(dpi = 900)

# Define geographic bounds for contiguous U.S.
xmin, xmax = -125, -65  
ymin, ymax = 24, 50      

# Function to generate each map (returns a plot object)
function plot_variable(data, variable_name, title_name)

    p = plot(color=:blues, grid=false, framestyle=:none, axis=nothing, ratio=1, 
             xlims=(xmin, xmax), ylims=(ymin, ymax), title=title_name, size = (900,675), colorbar_size=0.005)

            
min_val = minimum(skipmissing(data[!,variable_name]))
max_val = maximum(skipmissing(data[!,variable_name]))

# Define the custom 5-color scale
custom_cmap = cgrad([:red, :yellow, :green],  scale=:symlog) 

    for row in eachrow(data)
        if !ismissing(row[variable_name])
            plot!(p, row.geometry, fill_z=row[variable_name], label="", lw=0.5, linecolor=:black, c=custom_cmap)
        end
    end

    return p  # Ensure function returns a plot object
end

# Generate maps for the four variables
p1 = plot_variable(merged_df, :old_w, "Old Wages by County")
p2 = plot_variable(merged_df, :cf_hatL, "Change in Workers by County")
p3 = plot_variable(merged_df, :cf_hatR, "Change in Residents by County")
p4 = plot_variable(merged_df, :diag_old_λ, "Change in Living Where You Work by County")


# Save the plots to the output/figures folder with higher resolution
savefig(p1, "output/figures/old_wages_by_county.png")
savefig(p2, "output/figures/change_in_workers_by_county.png")
savefig(p3, "output/figures/change_in_residents_by_county.png")
savefig(p4, "output/figures/change_in_living_where_you_work_by_county.png")

