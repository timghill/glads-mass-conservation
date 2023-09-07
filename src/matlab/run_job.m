% function pa = run_job()
set_paths;
para_steady = get_para_steady();
pa = para_steady;
pa.physical
run_model(para_steady);