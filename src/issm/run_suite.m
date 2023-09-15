job_arr = { {"seasonal_001.mat", 0.0095, 5./4., 3./2., 0},
            {"seasonal_002.mat", 0.0189, 3./2., 3./2., 0},
            {"seasonal_003.mat", 0.1,  3.0,   2.0,   0},
            {"seasonal_004.mat", 0.1, 5./4., 3./2.,  1/2000},
            {"seasonal_005.mat", 0.1, 3./2., 3./2.,  1/2000}};

job_ids = [2, 3, 4, 5];
% job_ids = [1];
for jid=job_ids
    job_paras = job_arr{jid};
    fname = job_paras{1};
    k_s = job_paras{2};
    alpha_s = job_paras{3}; beta_s = job_paras{4};
    omega = job_paras{5};
    run_single_job(fname, k_s, alpha_s, beta_s, omega);
end
