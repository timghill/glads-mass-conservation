function pa = run_job(fname, k_s, alpha, beta, omega, k_c)
set_paths;
para_steady = get_para_seasonal(fname, k_s, alpha, beta, omega, k_c);
pa = para_steady;
pa.physical
run_model(para_steady);
end