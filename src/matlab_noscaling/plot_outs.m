outs = load('RUN/seasonal_001.mat');

dmesh = load('../data/mesh.mat');

tindex = 170;
figure('Units', 'inches', 'Position', [3, 2, 8, 3])
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes, ...
    'FaceVertexCData', outs.phis(:, tindex), 'EdgeColor', 'none', ...
    'FaceColor', 'flat');
axis image
cb = colorbar;
cb.Label.String = '\phi (MPa)';
print('phi_statevariables', '-dpng', '-r600')
% 
% figure
% plot(mean(outs.phis))
% cmocean('dense')