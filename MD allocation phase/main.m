format compact
clear
clc
clf reset

% ----------

instance = 1;

res = 8;
TA_start_range = linspace(0,1,res+1);
TA_start_range(end) = [];
omega_end_range = linspace(0,360,res+1);
omega_end_range(end) = [];

instance_sim_indexes = [35:53];

filename = "MD_allocation_instance_"+string(instance)+".mat";

if ~exist(filename, 'file')
    ind_s = 1;
    pso_results_struct = struct();
    save(filename,"ind_s","pso_results_struct")
end

max_iters = 4e2;
plot_result = true;

target_altitude = 1371.2e3;

load(filename);

while  ind_s <= numel(instance_sim_indexes)

    load(filename);

    [ind_o,ind_t] = ind2sub([numel(TA_start_range),numel(omega_end_range)],instance_sim_indexes(ind_s));

    TA_start_norm = TA_start_range(ind_t);
    target_omega = omega_end_range(ind_o);

    %dv_to_beat = -250*cosd(target_omega)+3.75e3;
    dv_to_beat = 3.75e3;

    if rem(target_omega,360) == 0
        dv_to_beat = 3e3;
    end

    PSO_failed = true;
    min_manuever_cycle = 0;
    failed_attempts = 0;

    while PSO_failed

        dv_kill_markers = [
        35, 5.5e3
        100, dv_to_beat + 750
        150, dv_to_beat + 500
        max_iters, dv_to_beat
        ]; %if not below right dv at left iter then the PSO session is killed

        [global_best_struct, PSO_failed, steps_forward] = PSO_handler(target_altitude, TA_start_norm, target_omega, max_iters, dv_kill_markers, plot_result, min_manuever_cycle);

        failed_attempts = failed_attempts+1;

        if failed_attempts > 0 && rem(target_omega,360) ~= 0 && steps_forward ~= 0
            dv_to_beat = dv_to_beat+300; %alas
            dv_to_beat = min([dv_to_beat,4.3e3]);
        end

        save(filename,"ind_s","pso_results_struct","target_altitude","failed_attempts") 
    end

    %save results
    pso_results_struct = setfield(pso_results_struct,"ind_s_"+string(instance_sim_indexes(ind_s)),"results_struct", global_best_struct);
    pso_results_struct = setfield(pso_results_struct,"ind_s_"+string(instance_sim_indexes(ind_s)),"TA_omega", [TA_start_norm*360, target_omega]);
    fprintf("completed sim %i.\n",ind_s);

    ind_s = ind_s+1;

    save(filename,"ind_s","pso_results_struct","target_altitude","failed_attempts") 
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);
