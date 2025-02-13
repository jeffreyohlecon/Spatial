function diff = fn_mode_111513(x,dist_rail,dist_road,dist_water,dist_air,frac_rail,frac_water,frac_air)

    % here are the estimates (all scaled by theta)
    cost_road_var = (x(1));
    cost_rail_var = (x(2));
    cost_water_var = (x(3));
    cost_air_var = (x(4));    
    cost_rail_fix = (x(5));
    cost_water_fix = (x(6));
    cost_air_fix = (x(7));

% calculating the method of moments 
       denom = (exp(-cost_road_var*dist_road))+...
           (exp(-cost_rail_var*dist_rail-cost_rail_fix))+...
           (exp(-cost_water_var*dist_water-cost_water_fix))+...
           (exp(-cost_air_var*dist_air-cost_air_fix));
       
       pred_frac_rail = (exp(-cost_rail_var*dist_rail-cost_rail_fix)) ./  denom;
       pred_frac_water = (exp(-cost_water_var*dist_water-cost_water_fix)) ./ denom;
       pred_frac_air = (exp(-cost_air_var*dist_air-cost_air_fix)) ./ denom;
       
% How close are we?
    diff = (nanmean(frac_rail - pred_frac_rail).^2).^0.5 + ...
        (nanmean(frac_water - pred_frac_water).^2).^0.5 + ...
        (nanmean(frac_air - pred_frac_air).^2).^0.5;
