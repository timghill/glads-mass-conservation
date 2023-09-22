function para = get_para(fname, k_s, alpha, beta, omega, k_c)
% para = get_para_steady(mesh_nr)
%
% Get default parameters for SHMIP ice-sheet margin domain run

addpath(genpath('../data/'))

para = get_default_para();  % Start with model defaults
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Model description
pm.model_run_descript = 'Run to steady';
pm.git_revision_model_runs = strtrim(git('rev-parse --verify HEAD'));
pm.verbosity = 10;   % Lots of output details

%% Model output directories
pm.dir.model_save = ['./', 'RUN', '/'];
pm.save_filename = [pm.dir.model_save, fname];
pm.save_filename_root = '';
pm.save_index_file = 0;

%% Mesh
pm.file.mesh = 'data/mesh.mat';
dmesh = load(pm.file.mesh);

%%  Physical parameters
pp.l_bed = 10;
pp.h_bed = 0.5;
pp.l_c = pp.l_bed;

pp.cond_s = k_s;
pp.alpha_s = alpha;
pp.beta_s = beta;
pp.omega = omega;

e_v = 1e-4;
pin.e_v = make_anon_fn('@(xy) double(0*xy(:,1) + e_v)',e_v);

pp.cond_c = k_c;

const_u_bed = 30/pp.year;
pin.u_bed = make_anon_fn('@(xy) double(0*xy(:, 1) + const_u_bed)',const_u_bed);

pp.creep_const_c = 1.78e-25;
pp.creep_const_s = 1.78e-25;
% pp.creep_const_s_soft = config.creep_const_soft; % Applied when N<0
pp.creep_const_s_soft = pp.creep_const_s; % Applied when N<0
% pp.creep_const_s_soft = 0;

%% Matlab-specific
pp.flags.max_S = 500;

op = {'freeze-on', 'freeze-melt', 'none'};
pp.flags.include_sheet_diss_press = op{1}; % if true, allows freezing shut of the sheet on reverse slopes


pp.float_frac = 0; % used below for BC

%% Numerics
pn.zero_channels_on_boundary = 1;
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{5};

st = {'ode15s', 'ode23s', 'ode23t', 'odebim'};  % can also use ode23t, ode23s but ode15s is probbaly best
pn.ts.ode15s.stepper = st{1};

% pn.ts.ode15s.opts = odeset('AbsTol', 1e-1, 'RelTol', 1e-2);

%% Domain geometry
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed(xy));');
pin.ice_thickness = make_anon_fn('@(xy, time) double(thick(xy));');

%% INPUT FUNCTIONS
% Boundary conditions

% Dirichlet BC for phi: applied at nodes where bmark is odd
pin.bc_dirichlet_phi = make_anon_fn('@(xy, time, bmark, phi_0, phi_m) double(pp.float_frac * (phi_0-phi_m) + phi_m)', pp);

% Flux BC for phi and h_w: i.e. set phi or h_w such that this flux is
% given. Applied at edges where bmark_edge is even
% zero flux:
pin.bc_flux = make_anon_fn('@(xy, time, bmark_edge) double(zeros(sum(~logical(mod(bmark_edge,2)) & bmark_edge>0),1))');

% Initial conditions
pin.ic_h = make_anon_fn('@(xy, time) double(0.1*pp.h_bed + 0*xy(:, 1));', pp);
pin.ic_S =  make_anon_fn('@(xy, time) double(0*xy(:,1));');

%% Partially nondimensionalize and wrap

% Set scales for physical constants
ps.rho_w = pp.rho_w;
ps.rho_i = pp.rho_i;
ps.L_fusion = pp.L_fusion;
ps.g_grav = pp.g_grav;

% ps.rho_w = 1;
% ps.rho_i = 1;
% ps.L_fusion = 1;
% ps.g_grav = 1;

% No scaling for dependent variables
ps.phi = 1;
ps.Q = 1;
ps.u_bed = 1;

% No scaling for spatial coordinates
ps.x = 1;
ps.x_offset = 0; ps.y_offset = 0;

% No scaling for model parameters
ps.cond_s = 1;
ps.cond_c = 1;
ps.l_bed = 1;
ps.h_bed = 1;
ps.l_c = 1;
ps.creep_const_s = 1;

ps.creep_const_c = ps.creep_const_s;
ps.h = ps.h_bed;
ps.h_w = ps.h;

% Compute scales for dependent variables (copied from set_default_scales)
if pp.omega>0
    % If we are using transition flux parameterization, pp.alpha_s
    % is interpreted as the flow exponent in the turbulent limit, but
    % we want to use the flow in the laminar limit to set the scale
    % for q
    alpha_scale = 3;
    beta_scale = 2;
else
    % We are using standard flux parameterization, use alpha_s and beta_s
    % as is
    alpha_scale = pp.alpha_s;
    beta_scale = pp.beta_s;
end
ps.q = ps.cond_s.*ps.h.^alpha_scale.*(ps.phi/ps.x).^(beta_scale-1); % sheet discharge (m^2/s) [q]
ps.nu = ps.q;       % kinematic viscosity of water (m2 s-1)
ps.l_spacing = ps.Q/ps.q; % average channel spacing (m) [l]
ps.S = ps.h_bed*ps.l_spacing; % channel cross sectional area (m^2) [S]
ps.S_w = ps.S; % channel water filled cross sectional area (m^2) [S_w]
ps.z = ps.phi/ps.rho_w/ps.g_grav; % vertical length scale for any elevation (m) [B] [H]
ps.t = ps.x*ps.h/ps.q; % time scale (s) [t]

%% Input fields scales
ps.source_term_s = ps.q/ps.x;  % source term sheet (ms^-1) [m]
ps.source_term_c = ps.Q;  % source term channel (m^3 s^-1) (this happens at the nodes) [Q_s or/and Q_m]
ps.e_v = ps.rho_w*ps.g_grav*ps.h/ps.phi; % englacial-void-ratio scale
ps.moulin_area = 1/( ps.phi/ps.rho_w/ps.g_grav/(ps.S*ps.x) ) ; % moulin area scale [A_m]

[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);

para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
