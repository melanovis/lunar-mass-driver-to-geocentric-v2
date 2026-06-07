function [jacobian, residual_matrix, r_target, state_earth_rendezvous] = fetch_jacobian_and_residual(r_target_offset, ER3BP_timespan, ER3BP_state_inital, lunar_mass, earth_mass, error_tolerance)

    peturb_factor = 5e-6;

    t_start = ER3BP_state_inital(1);

    [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow_harsh));
    scene_timerange = scene_timerange + t_start;
    scene_rocket_statematrix = [scene_timerange, scene_statematrix(:,1:6)];

    state_earth_rendezvous = scene_statematrix(end,7:12);
    r_target = scene_statematrix(end,7:9) + r_target_offset; %targeting some offset from the earth

    r_original = scene_rocket_statematrix(end,2:4);
    residual_matrix = (r_original - r_target).';

    jacobian = zeros(3,4);

    for n=1:4
        
        if n~=4
            ER3BP_state_test = ER3BP_state_inital;
            ER3BP_time_test = ER3BP_timespan;
            ER3BP_state_test(4+n) = ER3BP_state_test(4+n) + peturb_factor;
        else
            ER3BP_state_test = ER3BP_state_inital;
            ER3BP_time_test = ER3BP_timespan;
            ER3BP_time_test(2) = ER3BP_time_test(2) + peturb_factor;
        end

        [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_time_test, ER3BP_state_test ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow_harsh));
        scene_timerange = scene_timerange + t_start;
        scene_rocket_statematrix = [scene_timerange, scene_statematrix(:,1:6)];
        r_approach = scene_rocket_statematrix(end,2:4);
        jacobian(:,n) = ((r_approach - r_original)/peturb_factor).';
    end

end

function status = stop_too_slow_harsh(t,y,flag)
    persistent start_time
    status = 0;
    if strcmp(flag,'init')
        start_time = tic;
    elseif isempty(flag)
        time_max = 0.05; %time limit for cutting off ode solver
        if toc(start_time) > time_max
            status = 1; %stop ode solver
        end
    end
end

