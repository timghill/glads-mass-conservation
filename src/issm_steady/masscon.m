% Check mass conservation for ISSM model outputs
set_paths;
outs = load('MassCon.mat');
md = outs.md;
clear outs
% Add fields to md.mesh
% md = meshconvert(md,md.mesh.elements,md.mesh.x,md.mesh.y);

Q = [md.results.TransientSolution.ChannelDischarge];
time = [md.results.TransientSolution.time];
N = [md.results.TransientSolution.EffectivePressure];
S = [md.results.TransientSolution.ChannelArea];

% Externally load mesh
dmesh = load('../src/data/mesh.mat');

Q_node = sparse(dmesh.tri.n_nodes, dmesh.tri.n_edges);
for ii=1:dmesh.tri.n_nodes
    neigh_edges = dmesh.tri.connect_edge_inv{ii};
    num_neigh = length(neigh_edges);
    Q_signs = zeros(num_neigh, 1);
    for jj=1:num_neigh
        n1 = dmesh.tri.connect_edge(neigh_edges(jj), 1);
        n2 = dmesh.tri.connect_edge(neigh_edges(jj), 2);
        if n1==ii
            Q_signs(jj) = 1;
        elseif n2==ii
            Q_signs(jj) = -1;
        end
    end
    Q_neigh = Q(neigh_edges, end);
    Q_node(ii, neigh_edges) = Q_neigh.*Q_signs;
end

Q_sum = full(sum(Q_node, 2));
Q_sum_interior = Q_sum(dmesh.tri.bmark==0);
Q_sum_outlet = Q_sum(dmesh.tri.bmark==1);
n_nz = length(find(abs(Q_sum)>1e-3));

figure
hold on
scatter(1:length(Q_sum_interior), Q_sum_interior, 'Marker', '.')

figure
hold on
scatter(1:length(Q_sum_outlet), Q_sum_outlet, 'Marker', '.')

figure('Units', 'inches', 'Position', [2, 2, 8, 4])
T = tiledlayout(1, 1);
ax1 = nexttile;
hold on
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', Q_sum, 'FaceColor', 'flat')
cmocean('balance')
clim([-1, 1])
axis image
xlim([0, 100e3])
ylim([0, 25e3])
cb = colorbar;
cb.Label.String = 'Q_{\rm{node}} (m^3 s^{-1})';
title('ISSM-GlaDS')

Q_cmap = cmocean('turbid');
Qmin = 1e-1; Qmax = 100;
plotchannels(md, abs(md.results.TransientSolution(end).ChannelDischarge),...
        'colormap', Q_cmap, 'min', Qmin, 'max', Qmax)
axis image
xlim([0, 100e3])
ylim([0, 25e3])

% Add secondary invisible axes for the second colorbar
cax = axes(T);
cax.Visible=false;
colormap(cax, Q_cmap);
cb2 = colorbar(cax);
cb2.Layout.Tile = 'North';
cb2.Label.String = 'Q (m^3 s^{-1})';
cb2.Ticks = [1, 5, 10, 15, 20];
clim(cax, [Qmin, Qmax])

% print('ISSM_glads_q_net', '-dpng', '-r600')

figure
hold on
Qplot = abs(md.results.TransientSolution(end).ChannelDischarge);
Qplot(Qplot<Qmin) = nan;
edge_plot(gca, dmesh, Qplot, cmocean('turbid'), [Qmin, Qmax], 'vmin', 0,...
    'Lineweights', [2, 2])
axis image
xlim([0, 100e3])
ylim([0, 25e3])