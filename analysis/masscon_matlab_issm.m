% Compare mass conservation for matlab and ISSM GlaDS implementations

mat_fname = '../src/matlab/RUN/output.mat';
issm_fname = '../src/issm/MassCon_laminar.mat';
mesh_fname = '../src/data/mesh.mat';

dmesh = load(mesh_fname);
addpath('../src/data/')
[bed, thick] = bed_thickness(dmesh.tri.nodes);
phi_bed = 1e3*9.8*bed;
p_ice = 910*9.8*thick;

%% MATLAB
out = load(mat_fname);
out_pp = pp_do_pp(out);
Q_mat = out_pp.Q(:, end);

N_mat = out_pp.N(:, end);
phi_mat = out_pp.phi(:, end);
pw_mat = phi_mat - phi_bed;

N2_mat = p_ice - pw_mat;
ff_mat = pw_mat./p_ice;

%% ISSM
out_issm = load(issm_fname);
md = out_issm.md;
Q_issm = out_issm.md.results.TransientSolution(end).ChannelDischarge;
N_issm = out_issm.md.results.TransientSolution(end).EffectivePressure;
phi_issm = out_issm.md.results.TransientSolution(end).HydraulicPotential;
all_Q_issm = [md.results.TransientSolution.ChannelDischarge];
pw_issm = phi_issm - phi_bed;
N2_issm = p_ice - pw_issm;
ff_issm = pw_issm./p_ice;

Q_node_mat = sparse(dmesh.tri.n_nodes, dmesh.tri.n_edges);
Q_node_issm = sparse(dmesh.tri.n_nodes, dmesh.tri.n_edges);

for ii=1:dmesh.tri.n_nodes
    neigh_edges = dmesh.tri.connect_edge_inv{ii};
    for jj=1:length(neigh_edges)
        edgenum = neigh_edges(jj);
        neigh_nodes = dmesh.tri.connect_edge(edgenum, :);
        if neigh_nodes(1)==ii
            q_sign = 1;
        else
            q_sign = -1;
        end
        Q_node_mat(ii, edgenum) = Q_mat(edgenum).*q_sign;
        Q_node_issm(ii, edgenum) = Q_issm(edgenum).*q_sign;
    end
end

Q_net_mat = full(sum(Q_node_mat, 2));
Q_net_issm = full(sum(Q_node_issm, 2));

% Flotation
figure
subplot(2, 1, 1)
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
    'FaceVertexCData', ff_mat, 'EdgeColor', 'none', 'FaceColor', 'flat')
cmocean('haline')
clim([0, 1])
axis image
colorbar
title('Matlab')

subplot(2, 1, 2)
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
    'FaceVertexCData', ff_issm, 'EdgeColor', 'none', 'FaceColor', 'flat')
cmocean('haline')
clim([0, 1])
axis image
colorbar
title('ISSM')

% N
figure
subplot(2, 1, 1)
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
    'FaceVertexCData', N_mat, 'EdgeColor', 'none', 'FaceColor', 'flat')
cmocean('thermal')
clim([0, 2.5e6])
axis image
colorbar
title('Matlab')

subplot(2, 1, 2)
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
    'FaceVertexCData', N_issm, 'EdgeColor', 'none', 'FaceColor', 'flat')
cmocean('thermal')
clim([0, 2.5e6])
axis image
colorbar
title('ISSM')


% Mass conservation

fig = figure('units', 'inches', 'position', [4, 4, 8, 6]);
T = tiledlayout(2, 1);

ax1 = nexttile(1);

patch(ax1, 'Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
    'FaceVertexCData', Q_net_mat, 'EdgeColor', 'none', 'FaceColor', 'flat');
cmocean('balance')
axis image
title('MATLAB')
cb1 = colorbar;
cb1.Label.String = 'Q_{net}';
clim([-1, 1])
% ax,dmesh,Cdata,cmap,caxis,varargin)
Q_mat_plot = abs(Q_mat);
Q_mat_plot(Q_mat_plot<1) = nan;

cax1 = axes(T);
edge_plot(gca, dmesh, Q_mat_plot, cmocean('turbid'), [1, 100], 'vmin', 1)
cax1.Visible = false;
cb2 = cax1.Colorbar;
cb2.Layout.Tile = 'North';
cb2.Label.String = 'Q (m^3 s^{-1})';

ax2 = nexttile(2);

patch(ax2, 'Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
    'FaceVertexCData', Q_net_issm, 'EdgeColor', 'none', 'FaceColor', 'flat');
cmocean('balance')
% clim([-1, 1])
axis image
title('ISSM')
colorbar
clim([-1, 1])
Q_issm_plot = abs(Q_issm);
Q_issm_plot(Q_issm_plot<1) = nan;

cax2 = axes(T);
cax2.Layout.Tile = 2;
edge_plot(gca, dmesh, Q_issm_plot, cmocean('turbid'), [1, 100], 'vmin', 1)
colorbar(cax2, 'off')
cax2.Visible = false;

% print(fig, 'glads_node_conservation', '-dpng', '-r600')

% figure
% hold on
% plotchannels(md, abs(Q_issm), 'min', 1, 'max', 100, 'colormap', cmocean('turbid'))
% axis image
% xlim([0, 100e3])
% ylim([0, 25e3])
% 
% figure
% hold on
% edge_plot(gca, dmesh, abs(md.results.TransientSolution(end-30).ChannelArea), cmocean('turbid'), [1, 100], 'vmin', 1)

figure
hold on
scatter(N_mat, N_issm)
xlabel('MATLAB')
ylabel('ISSM')
axis image
xlim([0.5e6, 2.5e6])
ylim([0.5e6, 2.5e6])
grid on
plot([0, 1], [0, 1], 'k')
