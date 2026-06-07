format compact
clear
clc
clf reset

% ----------

instances = [1:7];

TA_start_range = linspace(0,1,8+1);
TA_start_range(end) = [];
omega_end_range = linspace(0,360,8+1);
omega_end_range(end) = [];

pso_results_total = struct();
for n=1:numel(instances)

    str_filename = "MD_allocation_instance_"+instances(n)+".mat";
    load(str_filename)

    ind_sim_fields = fieldnames(pso_results_struct);
    for ind_field = 1:numel(ind_sim_fields)
        fieldname_spec = string(ind_sim_fields(ind_field));
        pso_results_total.(fieldname_spec) = pso_results_struct.(fieldname_spec);

        TA_omega = pso_results_total.(fieldname_spec).TA_omega;
        ind_t = find(TA_omega(1)==TA_start_range.*360);
        ind_o = find(TA_omega(2)==omega_end_range);
        dv_matrix(ind_t,ind_o) = pso_results_total.(fieldname_spec).results_struct.dv_total;
        acc_series = pso_results_total.(fieldname_spec).results_struct.acc_series;
        TW_matrix(ind_t,ind_o) = max(acc_series)/9.81;
    end
end

save("pso_results_collected.mat","pso_results_total")