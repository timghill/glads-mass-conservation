job_arr = { {"seasonal_001.mat", 0.0095, 5./4., 3./2., 0,      0.1},
            {"seasonal_002.mat", 0.0189, 3./2., 3./2., 0,      0.1},
            {"seasonal_003.mat", 0.1,    3.0,   2.0,   0,      0.1},
            {"seasonal_004.mat", 0.1,    5./4., 3./2., 1/2000, 0.1},
            {"seasonal_005.mat", 0.1,    3./2., 3./2., 1/2000, 0.1},
            {"seasonal_101.mat", 0.0095, 5./4., 3./2., 0,      0.5},
            {"seasonal_102.mat", 0.0189, 3./2., 3./2., 0,      0.5},
            {"seasonal_103.mat", 0.1,    3.0,   2.0,   0,      0.5},
            {"seasonal_104.mat", 0.1,    5./4., 3./2., 1/2000, 0.5},
            {"seasonal_105.mat", 0.1,    3./2., 3./2., 1/2000, 0.5}};

% job_ids = [1, 2, 3, 4, 5];
job_ids = 7:9;
for jid=job_ids
    job_paras = job_arr{jid};
    fname = job_paras{1};
    k_s = job_paras{2};
    alpha_s = job_paras{3}; beta_s = job_paras{4};
    omega = job_paras{5};
    k_c = job_paras{6};
    run_single_job(fname, k_s, alpha_s, beta_s, omega, k_c);
end
