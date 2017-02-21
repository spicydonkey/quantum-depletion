
LIST_N_PER_SHOT=10.^[0:4];
for iter=1:numel(LIST_N_PER_SHOT)
    N_PER_SHOT=LIST_N_PER_SHOT(iter);
    run('k_dist_tester.m');
end