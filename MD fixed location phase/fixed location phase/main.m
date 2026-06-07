format compact
clear
clc
clf reset

% ----------

instance = 6;
instance_sim_indexes = [96:100];

record_video = false;

TA_start_range = linspace(0,1,27+1);
TA_start_range(end) = [];
omega_end_range = linspace(0,360,10+1);
omega_end_range(end) = [];

filename = "MD_fixed_instance_"+string(instance)+".mat";

if ~exist(filename, 'file')
    ind_s = 1;
    pso_results_struct = struct();
    save(filename,"ind_s","pso_results_struct")
end

max_iters = 2.5e2;
plot_result = true;

target_altitude = 1371.2e3;

view_struct = struct();
view_struct = setfield(view_struct,"view_az",45);
view_struct = setfield(view_struct,"view_el",20);

load(filename);

if record_video
    v = VideoWriter("record_opt_s"+string(ind_s), 'MPEG-4');
    v.FrameRate = 60;
    open(v);
else
    v = [];
end
recorder_struct = struct();
recorder_struct = setfield(recorder_struct,"video_obj",v);
recorder_struct = setfield(recorder_struct,"record",record_video);
recorder_struct = setfield(recorder_struct,"frame_width",1920*1.5);

while  ind_s <= numel(instance_sim_indexes)

    load(filename);

    [ind_o,ind_t] = ind2sub([numel(omega_end_range),numel(TA_start_range)],instance_sim_indexes(ind_s));

    TA_start_norm = TA_start_range(ind_t);
    target_omega = omega_end_range(ind_o);

    dv_to_beat = 3.5e3;

    PSO_failed = true;
    failed_attempts = 0;

    while PSO_failed

        dv_kill_markers = [
        40, 5e3
        100, 4.2e3
        max_iters, dv_to_beat
        ];

        [global_best_struct, PSO_failed, steps_forward, view_struct] = PSO_handler(target_altitude, TA_start_norm, target_omega, max_iters, dv_kill_markers, plot_result, recorder_struct, view_struct);

        failed_attempts = failed_attempts+1;

        if failed_attempts > 0 && steps_forward ~= 0
            dv_to_beat = dv_to_beat+300;
            dv_to_beat = min([dv_to_beat,4e3]);
        end

        save(filename,"ind_s","pso_results_struct","target_altitude","failed_attempts","view_struct") 
    end

    %save results
    pso_results_struct = setfield(pso_results_struct,"ind_s_"+string(instance_sim_indexes(ind_s)),"results_struct", global_best_struct);
    pso_results_struct = setfield(pso_results_struct,"ind_s_"+string(instance_sim_indexes(ind_s)),"TA_omega", [TA_start_norm*360, target_omega]);
    fprintf("completed sim %i.\n",ind_s);

    ind_s = ind_s+1;

    save(filename,"ind_s","pso_results_struct","target_altitude","failed_attempts") 
end

if record_video
    close(v);
end

sound(sin(2*pi*400*(0:1/14400:0.15)), 14400);