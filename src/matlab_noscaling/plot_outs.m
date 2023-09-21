outs = load('RUN/seasonal_001_norho.mat');

dmesh = load('../data/mesh.mat');

tindex = 170;
f = figure('Units', 'inches', 'Position', [3, 2, 8, 3]);
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes, ...
    'FaceVertexCData', outs.phis(:, tindex), 'EdgeColor', 'none', ...
    'FaceColor', 'flat');
axis image
cb = colorbar;
cb.Label.String = '\phi (MPa)';
print(f, 'phi_norho', '-dpng', '-r600')

f2 = figure;
plot(mean(outs.phis))
cmocean('dense')
print(f2, 'phi_norho_time', '-dpng', '-r600')