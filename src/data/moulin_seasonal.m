function m_moulin = moulin_inputs(meshtri, num_moulins, t)
    % Moulin inputs (m3.s-1) for GlaDS simulations. t measured in YEARS,
    % moulin input rate is m3.s-1

    % CONSTANTS
    amp_input = 5;
    mean_input = 0;

    % Randomly choose num_moulins interior nodes
    interior = find(meshtri.bmark==0);
    rng(40);
    interior_indices = randi(length(interior), num_moulins, 1);
    moulin_indices = interior(interior_indices);
    
    % Compute moulin inputs
    m_moulin = zeros(meshtri.n_nodes, length(t));
    seasonal = max(mean_input - cos(2*pi*t), 0);
    m_moulin(moulin_indices, :) = amp_input.*repmat(seasonal, [num_moulins,1]);
end
