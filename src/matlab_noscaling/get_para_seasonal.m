function para = get_para_steady(fname, k_s, alpha, beta, omega, k_c)
% para = get_para_steady(config)
%
% Set para for steady state run

%% Get defaults and unwrap
addpath('../')
para = get_para(fname, k_s, alpha, beta, omega, k_c);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Time
pt.end   = 5*pp.year;  % end time
pt.out_t = pt.start:10*86400:pt.end;

addpath(genpath('../data'))
pin.source_term_s = make_anon_fn('@(xy, time) double(0.05/86400/365 + 0*xy(:, 1));');
pin.source_term_c = make_anon_fn('@(time) double(moulin_seasonal(dmesh.tri, 50, time./pp.year));', dmesh, pp);

%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
