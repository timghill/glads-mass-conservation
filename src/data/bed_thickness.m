function [bed, thickness] = bed_thickness(xy)
    % Compute bed elevation and ice thickness for synthetic topo
    bed = 350 + 0*xy(:, 1);
    thickness = 50 + 6*(sqrt(xy(:, 1) + 5e3) - sqrt(5e3));
end