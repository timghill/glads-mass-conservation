function t = thickness(xy)
    % Compute bed elevation and ice thickness for synthetic topo
    t = 50 + 6*(sqrt(xy(:, 1) + 5e3) - sqrt(5e3));
end