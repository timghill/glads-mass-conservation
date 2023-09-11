function m_moulin = moulin_inputs(meshtri, num_moulins, t)
    % Moulin inputs (m3.s-1) for GlaDS simulations
    interior = find(meshtri.bmark==0);
    rng(40);
    interior_indices = randi(length(interior), num_moulins, 1);
    moulin_indices = interior(interior_indices);
    m_moulin = zeros(meshtri.n_nodes, 1);
    m_moulin(moulin_indices, :) = 5;
end