function para = get_para()
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
pm.save_filename = [pm.dir.model_save, 'output.mat'];
pm.save_filename_root = '';
pm.save_index_file = 0;

%% Mesh
pm.file.mesh = 'data/mesh.mat';
dmesh = load(pm.file.mesh);

%%  Physical parameters
pp.cond_s = 0.01;
% pp.omega = 1/2000;

pp.l_bed = 2;
pp.h_bed = 0.1;

% pp.alpha_s = 3/2;
% pp.beta_s = 2;
pp.alpha_s = 5/4;
pp.beta_s = 3/2;
pp.omega = 0;

e_v = 1e-4;
pin.e_v = make_anon_fn('@(xy) double(0*xy(:,1) + e_v)',e_v);

pp.cond_c = 0.1;

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
pp.flags.include_sheet_diss_press = op{3}; % if true, allows freezing shut of the sheet on reverse slopes


pp.float_frac = 0; % used below for BC

%% Numerics
pn.zero_channels_on_boundary = 1;
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{5};

st = {'ode15s', 'ode23s', 'ode23t', 'odebim'};  % can also use ode23t, ode23s but ode15s is probbaly best
pn.ts.ode15s.stepper = st{1};

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

%% Nondimensionalize and wrap
ps = set_default_scales(ps, pp, dmesh);
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);

para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
