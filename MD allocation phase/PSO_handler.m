function [global_best_struct, PSO_killed,steps_forward] = PSO_handler(target_altitude, TA_start_norm, target_omega, max_iters, dv_kill_markers, plot_result, min_manuever_cycle)

    G = 6.6743e-11;
    error_tolerance = 1e-12;
    
    lunar_radius = 1.7374e6;
    earth_radius = 6.371e6; %m
    
    target_radius = target_altitude + earth_radius;
    
    lunar_obliquity = 6.68;
    earth_obliquity = 23.44;
        
    lunar_mass = 7.35e22; %kg
    earth_mass = 5.972e24;
    mu = G*(earth_mass + lunar_mass);
    
    approx_semimajor = 385692.5e3;
    
    pair_eccentricity = 0.0549;
    
    lunar_inc = 5.145;
    lunar_right_ascension = 0;
    lunar_argument_peri = 0;
    
    earth_inc = lunar_inc;
    earth_right_ascension = 0;
    earth_argument_peri = 0;
    
    %% defining earth-luna orbit
    
    %finding velocity magnitudes
    r_relative = approx_semimajor*(1-pair_eccentricity);
    v_relative = sqrt(mu*(1 + pair_eccentricity)/(approx_semimajor*(1 - pair_eccentricity)));
    
    r_earth_mag = (lunar_mass / (earth_mass + lunar_mass))*r_relative;
    r_luna_mag = (earth_mass / (lunar_mass + earth_mass))*r_relative;
    v_earth_mag = (lunar_mass / (earth_mass + lunar_mass))*v_relative;
    v_luna_mag = -(earth_mass / (earth_mass + lunar_mass))*v_relative;
    
    lunar_T_matrix = [
    cosd(lunar_right_ascension), -sind(lunar_right_ascension), 0
    sind(lunar_right_ascension), cosd(lunar_right_ascension), 0
    0, 0, 1
    ]*[
    1, 0, 0
    0, cosd(lunar_inc), -sind(lunar_inc)
    0, sind(lunar_inc), cosd(lunar_inc)
    ]*[
    cosd(lunar_argument_peri), -sind(lunar_argument_peri), 0
    sind(lunar_argument_peri), cosd(lunar_argument_peri), 0
    0, 0, 1
    ];
    
    lunar_r_initial = lunar_T_matrix*[0; r_luna_mag; 0];
    lunar_v_initial = lunar_T_matrix*[v_luna_mag; 0; 0];
    
    earth_T_matrix = [
    cosd(earth_right_ascension), -sind(earth_right_ascension), 0
    sind(earth_right_ascension), cosd(earth_right_ascension), 0
    0, 0, 1
    ]*[
    1, 0, 0
    0, cosd(earth_inc), -sind(earth_inc)
    0, sind(earth_inc), cosd(earth_inc)
    ]*[
    cosd(earth_argument_peri), -sind(earth_argument_peri), 0
    sind(earth_argument_peri), cosd(earth_argument_peri), 0
    0, 0, 1
    ];
    
    earth_r_initial = earth_T_matrix*[0; r_earth_mag; 0];
    earth_v_initial = earth_T_matrix*[v_earth_mag; 0; 0];
    
    earth_r_initial = -earth_r_initial;
    
    earthlunar_state_intial = [earth_r_initial.', earth_v_initial.', lunar_r_initial.', lunar_v_initial.'];
    
    earthluna_orbital_period = sqrt( (approx_semimajor^3)*(4*pi^2)/mu);
    
    timespan = linspace(0,earthluna_orbital_period*2,7e3);
    
    [timerange_init, statematrix_init] = ode113(@(timerange, statematrix) earth_moon_ODE(timerange, statematrix, earth_mass, lunar_mass), timespan, earthlunar_state_intial ,odeset('Reltol',error_tolerance));
    earth_init_statematrix = [timerange_init, statematrix_init(:,1:6)];
    lunar_init_statematrix = [timerange_init, statematrix_init(:,7:12)];
    
    lunar_peri_state = lunar_init_statematrix(1,:);
    earth_apo_state = earth_init_statematrix(1,:);
    
    lunar_TA_list = get_planet_TAlist(lunar_peri_state(2:4),lunar_peri_state(5:7),lunar_init_statematrix);
    earth_TA_list = get_planet_TAlist(earth_apo_state(2:4),earth_apo_state(5:7),earth_init_statematrix);
    earth_TA_list = rem(earth_TA_list+180,360); %TA function is designed to start at periapse, accounting for earth's apo start
    
    earth_TA_list = accumulate_increasing_TAlist(earth_TA_list);
    lunar_TA_list = accumulate_increasing_TAlist(lunar_TA_list);
    
    earth_init_statematrix = [earth_init_statematrix, earth_TA_list.'];
    lunar_init_statematrix = [lunar_init_statematrix, lunar_TA_list.'];
    
    
    %% PSO
    n_DOF = 21;
    var_size = [1, n_DOF]; %solution matrix size
    var_min = 0;
    var_max = 1;
    
    population = 8*round(interp1([0,1],[13,25],rand()));
    phi_1 = 2.05;
    phi_2 = 2.05;
    phi = phi_1+phi_2;
    kappa = 1;
    chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
    w = chi;
    w_damp = 0.75;
    w_original = 1;

    c1 = chi*phi_1;
    c2 = chi*phi_2;
    % c1 = 1;
    % c2 = 1;
    
    max_velocity = (var_max-var_min)*0.3;
    min_velocity = -max_velocity;
    
    %initalise particle template
    empty_particle.position = [];
    empty_particle.velocity = [];
    empty_particle.fitness = [];
    empty_particle.best.position = [];
    empty_particle.best.fitness = [];
    
    %initalise global best
    
    particle_init_series = [];
    particle_init_block = repmat(empty_particle, population, 1);
    particle = repmat(empty_particle, population, 1);
    
    if ~gcp().Connected 
        delete(gcp('nocreate'));
        parpool('local',8);
    end
    
    checkgood_threshold = round(interp1([0,1],[20,40],rand()));
    init_checkgood_quantity = 0;
    checkgood_series = [];
    
    while init_checkgood_quantity < checkgood_threshold
    
        checkgood_block = zeros(1,population);
        parfor n=1:population
            particle_init_block(n).position = unifrnd(var_min,var_max,var_size);
            input_vector = particle_init_block(n).position;
        
            [return_struct] = transfer_handler(input_vector, TA_start_norm, target_radius, target_omega, earth_init_statematrix, lunar_init_statematrix, earthluna_orbital_period, min_manuever_cycle);
        
            particle_init_block(n).fitness = return_struct.fitness;
            particle_init_block(n).velocity = zeros(var_size); 
            particle_init_block(n).best.position = particle_init_block(n).position;
            particle_init_block(n).best.fitness = particle_init_block(n).fitness;
            checkgood_block(n) = all(return_struct.check_matrix);
        end
    
        particle_init_series = [particle_init_series; particle_init_block];
        checkgood_series = [checkgood_series, checkgood_block];
        init_checkgood_quantity = sum(checkgood_series);
    end
    
    for n=1:length(particle_init_series(:))
        fit_init_list(n) = particle_init_series(n).fitness;
    end
    
    [~,ind_sortfit] = sort(fit_init_list);
    ind_sortfit = flip(ind_sortfit);
    
    ind_conv = find(checkgood_series==1);
    ind_sortfit(ind_conv) = [];
    ind_sortfit = [ind_conv, ind_sortfit];
    
    for n=1:population
        particle(n) = particle_init_series(ind_sortfit(n));
    end
    
    global_best.fitness = -inf;
    
    for n=1:population
        if particle(n).best.fitness > global_best.fitness
            global_best.position = particle(n).best.position;
            global_best.fitness = particle(n).best.fitness;
        end
        particle(n).velocity = rand(var_size).*max_velocity;
    end
    
    fprintf("------------------\n")
    
    iter = 1;
    ag_countdown = 10;
    global_best_fitness = -inf;
    global_best_dv = inf;
    global_best_dt = inf;
    steps_forward = 0;
    
    PSO_killed = false;
    
    global_best_struct = struct();
    
    while iter <= max_iters
    
        parfor n=1:population
            
            particle(n).velocity = w*particle(n).velocity ...
                + c1*rand(var_size).*(particle(n).best.position - particle(n).position) ...
                + c2*rand(var_size).*(global_best.position - particle(n).position);
            
            particle(n).velocity = max(particle(n).velocity, min_velocity);
            particle(n).velocity = min(particle(n).velocity, max_velocity);
        
            particle(n).position = particle(n).position + particle(n).velocity;
            
            particle(n).position = max(particle(n).position, var_min);
            particle(n).position = min(particle(n).position, var_max);
        
            input_vector = particle(n).position;
            [return_struct] = transfer_handler(input_vector, TA_start_norm, target_radius, target_omega, earth_init_statematrix, lunar_init_statematrix, earthluna_orbital_period, min_manuever_cycle);
        
            particle(n).fitness = return_struct.fitness;
        
            %update personal best
            if particle(n).fitness > particle(n).best.fitness
                particle(n).best.position = particle(n).position;
                particle(n).best.fitness = particle(n).fitness;
            end
           
        end
    
        % update_plot = false;
        % for n=1:population
        %     if particle(n).best.fitness > global_best_fitness
        %         input_vector = particle(n).best.position;
        %         [return_struct] = transfer_handler(input_vector, TA_start_norm, target_radius, target_omega, earth_init_statematrix, lunar_init_statematrix, earthluna_orbital_period, min_manuever_cycle);
        %         if all(return_struct.check_matrix) && return_struct.fitness > global_best_fitness
        %             %update global best
        %             global_best_struct = return_struct;
        %             global_best_fitness = global_best_struct.fitness;
        %             global_best_dv = global_best_struct.dv_total;
        %             global_best_dt = global_best_struct.dt_total;
        %             update_plot = true;
        %             steps_forward = steps_forward+1;
        %             global_best.position = particle(n).best.position;
        %             global_best.fitness = particle(n).best.fitness;
        %         end
        %     end
        % end
        
        update_plot = false;
        for n=1:population
            if particle(n).best.fitness > global_best.fitness
                global_best = particle(n).best;
                input_vector = particle(n).best.position;
                [return_struct] = transfer_handler(input_vector, TA_start_norm, target_radius, target_omega, earth_init_statematrix, lunar_init_statematrix, earthluna_orbital_period, min_manuever_cycle);
                if all(return_struct.check_matrix) && return_struct.fitness > global_best_fitness
                    %update global best
                    global_best_struct = return_struct;
                    global_best_fitness = global_best_struct.fitness;
                    global_best_dv = global_best_struct.dv_total;
                    global_best_dt = global_best_struct.dt_total;
                    update_plot = true;
                    steps_forward = steps_forward+1;
                end
            end
        end

        bestfitnesss(iter) = global_best.fitness;
        
        if iter > 20
            if std(bestfitnesss(end-3:end)) < 5e-3
                ag_countdown = ag_countdown-1;
            else
                ag_countdown = 10;
            end
            if ag_countdown <= 0
                fitness_score = [];
                for n=1:population
                    fitness_score(n) = particle(n).fitness;
                end
                [~,ind_sort] = sort(fitness_score);
        
                if rand() < 0.5
                    ind_ag = ind_sort(1:round(length(ind_sort)*0.5));
                else
                    ind_ag = ind_sort;
                end
    
                for n=1:length(ind_ag)
                    particle(ind_ag(n)).position = rand(1,n_DOF)*var_max;
                    particle(ind_ag(n)).velocity = rand(1,n_DOF)*max_velocity;
                end
        
                w = w_original;
                ag_countdown = 10;
                fprintf("\n agitating.\n")
            end
        end
        
        w = w * w_damp;
    
        fprintf("\n (%3.2f, %3.2f) completed iter %i, %3.4f, %3.6f km/s, %i.\n",TA_start_norm, target_omega, iter, global_best_fitness , global_best_dv / 1e3, steps_forward )
    
        if update_plot && plot_result

            cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3 ));

            scatter(nan,nan,"w")
            hold on
            scatter(nan,nan,"w")
            grid on
            axis tight equal
            plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b")
            plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"k")
    
            transit_statematrix = global_best_struct.transit_statematrix;
            patch([transit_statematrix(:,2);nan],[transit_statematrix(:,3);nan],[transit_statematrix(:,4);nan],[transit_statematrix(:,end);nan],'EdgeColor','interp',linewidth=1.5)
            colormap(cmap)
            clim([1,max(transit_statematrix(:,end))+1])
    

            %plot target orbit
            state_earth_arrival = transit_statematrix(end,8:13);
            state_luna_depart = transit_statematrix(1,14:19);
            
            target_orbit_stateinital = [global_best_struct.targeting_rv(1,:), global_best_struct.targeting_rv(2,:)] ;
            target_orbit_timespan = [0,sqrt((target_radius^3)*(4*pi^2)/(G*earth_mass))];
            [target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),target_orbit_timespan,target_orbit_stateinital,odeset('Reltol',error_tolerance));
            target_orbit_statematrix = [target_timerange,target_state_matrix];
            for n=2:4
                target_orbit_statematrix(:,n) = + target_orbit_statematrix(:,n) + state_earth_arrival(n-1);
            end
            [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(state_earth_arrival(1),state_earth_arrival(2),state_earth_arrival(3),earth_radius,earth_radius,earth_radius,64);
            earth_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[0,0,1]);
            plot3(target_orbit_statematrix(:,2),target_orbit_statematrix(:,3),target_orbit_statematrix(:,4),"k")

            [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(state_luna_depart(1),state_luna_depart(2),state_luna_depart(3),lunar_radius,lunar_radius,lunar_radius,64);
            luna_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[0.5,0.5,0.5]);

            legend("$\Delta$v: "+string(round(global_best_dv./1e3,3))+" km/s", "$\Delta$t: " + string(round(global_best_dt / (3600*24), 3 ))+" days" , Interpreter="latex", FontSize=20, location = "northwest")
            legend boxoff
    
            set(gcf, 'Color', [1,1,1])
            set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    
            % xlim([-4.5e8,4.5e8])
            % ylim([-4.5e8,4.5e8])
            % zlim([-2e8,2e8])

            view([45,20])
            drawnow()
            hold off
    
        end
    
        iter = iter+1;
    
        % PSO kill conditions
        if iter > 15 && isinf(global_best_fitness)
            PSO_killed = true;
            break
        end
    
        if ismember(iter,dv_kill_markers(:,1))
            ind_killmarker = find(iter==dv_kill_markers(:,1));
            if global_best_dv > dv_kill_markers(ind_killmarker,2)
                PSO_killed = true;
                break
            end
        end

    end
end



function TA_list_out = accumulate_increasing_TAlist(TA_list_in)
    a = 0;
    TA_prev = TA_list_in(1);
    for n=1:numel(TA_list_in)
        TA_curr = TA_list_in(n);
        if abs(TA_curr - TA_prev) > 180 %assume this means a rollback from 360 to 0
            a = a+360;
        end
        TA_list_out(n) = TA_curr + a;
        TA_prev = TA_curr;
    end
end

function TA_list = get_planet_TAlist(a,b,statematrix)
    a = a/norm(a);
    b = b/norm(b);

    plane_normal = cross(a,b);

    if plane_normal(3) < 0
        plane_normal = -plane_normal;
    end

    for n=1:height(statematrix)
        c = statematrix(n,2:4);
        TA_list(n) = rem( atan2d( dot(plane_normal,cross(a,c)) , dot(a,c) ) + 360, 360);
    end
end

function state_out = earth_moon_ODE(t, state_in, earth_mass, lunar_mass)

    G = 6.6743e-11;

    r_earth = state_in(1:3).';
    v_earth = state_in(4:6).';
    r_luna = state_in(7:9).';
    v_luna = state_in(10:12).';

    r_earthluna = r_earth - r_luna;
    r_unit_earthluna = norm(r_earthluna);

    acc_earth =  -G * lunar_mass * r_earthluna / r_unit_earthluna^3;
    acc_luna = G * earth_mass * r_earthluna / r_unit_earthluna^3;

    %state in derived for the next timestep
    state_out = [v_earth, acc_earth, v_luna, acc_luna].';
end

function state_out = earth_twobody_ODE(t,state,primary_mass)
    G = 6.67e-11;
    r = state(1:3); 
    v = state(4:6);
    r_mag = norm(r);
    r_unit = r/r_mag;
    force_gravity = r_unit*(-G*primary_mass/(r_mag^2));
    acc = force_gravity;
    state_out = [v; acc];
end
