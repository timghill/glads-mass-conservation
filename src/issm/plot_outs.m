set_paths;
outs = load('seasonal_005.mat');
md = outs.md;
clear outs

Q = [md.results.TransientSolution.ChannelDischarge];
time = [md.results.TransientSolution.time];
N = [md.results.TransientSolution.EffectivePressure];
S = [md.results.TransientSolution.ChannelArea];
% Externally load mesh
dmesh = load('../data/mesh.mat');
tindex = 164;


figure('Units', 'inches', 'Position', [2, 2, 8, 4])
T = tiledlayout(1, 1);
ax1 = nexttile;
hold on
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes, 'EdgeColor', 'none',...
        'FaceVertexCData', N(1:end, tindex), 'FaceColor', 'flat')
cmocean('balance')
clim([-3e6, 3e6])
axis image
xlim([0, 100e3])
ylim([0, 25e3])
cb = colorbar;
cb.Label.String = 'N (MPa)';
title('ISSM-GlaDS')

Q_cmap = cmocean('turbid');
Qmin = 1e-1; Qmax = 100;
plotchannels(md, abs(md.results.TransientSolution(tindex).ChannelDischarge),...
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