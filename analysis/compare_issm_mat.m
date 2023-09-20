cases = [1, 2, 3, 4, 5];
figname = 'issm_matlab';

% cases = [101, 102, 103, 104, 105];
% figname = 'issm_matlab_high_ks';

mat_fpattern = '../src/matlab_noscaling/RUN/seasonal_%03d.mat';
issm_fpattern = '../src/issm/seasonal_%03d.mat';
mesh_fname = '../src/data/mesh.mat';



dmesh = load(mesh_fname);
addpath('../src/data/')
[bed, thick] = bed_thickness(dmesh.tri.nodes);
phi_bed = 1e3*9.8*bed;
p_ice = 910*9.8*thick;

% Start some figures
fig_avg = figure;
ax_avg = axes(fig_avg);
hold on
grid on

node_mask = dmesh.tri.nodes(:, 1)>=27.5e3 & dmesh.tri.nodes(:, 1)<=32.5e3;
mat_linestyle = '-';
issm_linestyle = '--';
colors = [  0.579, 0.677, 0.781, 1;
           0.199, 0.328, 0.492, 1;
           0.250, 0.250, 0.250, 1;
           0.929, 0.835, 0.408, 1;
           0.836, 0.590, 0.160, 1];


for ii=1:length(cases)
    cid = cases(ii);
    mat_fname = sprintf(mat_fpattern, cid);
    issm_fname = sprintf(issm_fpattern, cid);
    
    %% MATLAB
    out = load(mat_fname);
    out_pp = pp_do_pp(out);
    Q_mat = out_pp.Q;
    
    N_mat = out_pp.N;
    phi_mat = out_pp.phi;
    pw_mat = phi_mat - phi_bed;
    
    N2_mat = p_ice - pw_mat;
    ff_mat = pw_mat./p_ice;

    tt_mat = out_pp.para.time.out_t/86400/365;
    
    %% ISSM
    out_issm = load(issm_fname);
    md = out_issm.md;
    Q_issm = [out_issm.md.results.TransientSolution.ChannelDischarge];
    N_issm = [out_issm.md.results.TransientSolution.EffectivePressure];
    phi_issm = [out_issm.md.results.TransientSolution.HydraulicPotential];
    pw_issm = phi_issm - phi_bed;
    N2_issm = p_ice - pw_issm;
    ff_issm = pw_issm./p_ice;
    tt_issm = [out_issm.md.results.TransientSolution.time];
    
    %% Plots
    
    % Spatial mean N
    
    figure
    hold on
    plot(mean(N_mat))
    plot(mean(N_issm))
    grid on
    legend({'Matlab', 'ISSM'})
    ylabel('N (MPa)')
    title(sprintf('Case %d', cid));
    
    % Snapshot flotation and channels
    tindex = 170;
    
    fig = figure('units', 'inches', 'position', [4, 4, 8, 6]);
    T = tiledlayout(2, 1);
    
    ax1 = nexttile(1);
    
    patch(ax1, 'Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
        'FaceVertexCData', ff_mat(:, tindex), 'EdgeColor', 'none', 'FaceColor', 'flat');
    cmocean('dense')
    axis image
    title(sprintf('MATLAB case %d', cid));
    cb1 = colorbar;
    cb1.Label.String = 'Q_{net}';
    clim([0, 1])
    % ax,dmesh,Cdata,cmap,caxis,varargin)
    Q_mat_plot = abs(Q_mat(:, tindex));
    Q_mat_plot(Q_mat_plot<1) = nan;
    
    cax1 = axes(T);
    edge_plot(gca, dmesh, Q_mat_plot, cmocean('turbid'), [1, 100], 'vmin', 1)
    cax1.Visible = false;
    cb2 = cax1.Colorbar;
    cb2.Layout.Tile = 'North';
    cb2.Label.String = 'Q (m^3 s^{-1})';
    
    ax2 = nexttile(2);
    
    patch(ax2, 'Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
        'FaceVertexCData', ff_issm(:, tindex), 'EdgeColor', 'none', 'FaceColor', 'flat');
    cmocean('dense')
    clim([0, 1])
    axis image
    title('ISSM')
    colorbar
    Q_issm_plot = abs(Q_issm(:, tindex));
    Q_issm_plot(Q_issm_plot<1) = nan;
    
    cax2 = axes(T);
    cax2.Layout.Tile = 2;
    edge_plot(gca, dmesh, Q_issm_plot, cmocean('turbid'), [1, 100], 'vmin', 1)
    colorbar(cax2, 'off')
    cax2.Visible = false;

    plot(ax_avg, tt_mat, mean(ff_mat(node_mask, :)), 'LineStyle', mat_linestyle,...
        'Color', colors(ii, :))
    plot(ax_avg, tt_issm, mean(ff_issm(node_mask, :)),'LineStyle', issm_linestyle,...
        'Color', colors(ii, :))
end

xlim(ax_avg, [4, 5])
xlabel(ax_avg, 'Year')
ylabel(ax_avg, 'p_w / p_i')

print(fig_avg, figname, '-dpng', '-r600')