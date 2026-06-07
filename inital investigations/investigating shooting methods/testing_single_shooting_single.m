format compact
clear
clc
clf reset

% ----------

G = 6.6743e-11;
error_tolerance = 1e-10;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; %m plus atmospheric margin

collision_tolerance = 100e3; %m, the border around objects to assume no collision

lunar_mass = 7.35e22; %kg
earth_mass = 5.972e24;
mu = G*(earth_mass + lunar_mass);

approx_semimajor = 385692.5e3;

pair_eccentricity = 0.0549;

lunar_SOI_radius = approx_semimajor*(1-pair_eccentricity)*(lunar_mass / (3*(earth_mass+lunar_mass)) ) ^ (1/3);

lunar_obliquity = 6.68;
earth_obliquity = 23.44;

lunar_inc = 5.145;
lunar_right_ascension = 0;
lunar_argument_peri = 0;

earth_inc = lunar_inc;
earth_right_ascension = 0;
earth_argument_peri = 0;

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

earth_init_statematrix = [earth_init_statematrix, earth_TA_list.'] ; %these are only for finding inital params
lunar_init_statematrix = [lunar_init_statematrix, lunar_TA_list.'] ;


%% building target orbit

target_radius = 1371.2e3 + earth_radius;
target_TA = 0;
target_inc = earth_obliquity + lunar_inc;
target_omega = 0; %controlled for precession

target_T_matrix = [
1, 0, 0
0, cosd(-earth_obliquity), -sind(-earth_obliquity)
0, sind(-earth_obliquity), cosd(-earth_obliquity)
]*[
cosd(target_omega), -sind(target_omega), 0
sind(target_omega), cosd(target_omega), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(target_inc), -sind(target_inc)
0, sind(target_inc), cosd(target_inc)
]*[
cosd(target_TA), -sind(target_TA), 0
sind(target_TA), cosd(target_TA), 0
0, 0, 1
];

target_r_initial = target_T_matrix*[0; target_radius; 0]; %this is what we're actually targeting when displaced to the earth statematrix
target_v_initial = target_T_matrix*[-sqrt(G*earth_mass/target_radius); 0; 0]; %in earth barycenter reference, not earth-moon barycenter reference

r_target_offset = target_r_initial.';
v_target_offset = target_v_initial.';

%for display only to make sure we haven't fucked it up
target_orbital_period = sqrt(((target_radius^3)/mu)*(2*pi)^2);
target_orbit_timespan = [0, target_orbital_period];
target_state_initial = [target_r_initial,target_v_initial];
[target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),target_orbit_timespan,target_state_initial,odeset('Reltol',error_tolerance));
target_orbit_statematrix = [target_timerange,target_state_matrix];

%% shooting method parameters

shooting_good = false;

while ~shooting_good

    r_shooting_start = [interp1([0,1],[-9,9],rand())*10^interp1([0,1],[6,7.5],rand()), interp1([0,1],[-9,9],rand())*10^interp1([0,1],[6,7.5],rand()), 0];
    
    ER3BP_timespan = [0,earthluna_orbital_period];
    
    %% regarding threebody sim
    
    TA_start_norm = 0;
    TA_start = rem(TA_start_norm*360,360);
    lunar_start_statematrix = interp_statematrix_via_TA(lunar_init_statematrix,TA_start);
    
    t_start = lunar_start_statematrix(1); %only for testing!
    
    earth_start_statematrix = interp_statematrix_via_time(earth_init_statematrix,t_start);
    
    [v_1,v_2,delta_t] = lambert_solve_anticlockwise(r_shooting_start, [0, -target_radius ,0], 0.1, G*earth_mass);
    v_guess_inital = v_1;
    
    rocket_start_statematrix = [t_start, r_shooting_start, v_guess_inital];
    ER3BP_state_inital = [rocket_start_statematrix(2:7), earth_start_statematrix(2:7), lunar_start_statematrix(2:7)];
    [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
    scene_timerange = scene_timerange + t_start;
    scene_rocket_statematrix = [scene_timerange, scene_statematrix(:,1:6)];
    
    ER3BP_timespan = t_start + [0, delta_t];
    
    hold on
    grid on
    
    shooting_good = false;
    
    tic
    for n=1:30
    
        [jacobian, residual_matrix, r_target, state_earth_rendezvous] = fetch_jacobian_and_residual(r_target_offset, ER3BP_timespan, ER3BP_state_inital, lunar_mass, earth_mass, error_tolerance);
        
        %fprintf("-------\n")
        norm(residual_matrix)
        
        [u_mat,s_mat,v_mat] = svd(jacobian);
        step_matrix = -v_mat * (s_mat \ (u_mat' * residual_matrix));
        v_step = step_matrix(1:3);
        t_step = step_matrix(4);
    
        ER3BP_state_inital(5:7) = ER3BP_state_inital(5:7) + v_step.';
        ER3BP_timespan(2) = ER3BP_timespan(2) + t_step;
    
        if norm(residual_matrix) < 5
            shooting_good = true;
            
            ER3BP_state_inital = single(ER3BP_state_inital);
            ER3BP_timespan = single(ER3BP_timespan);

            %just plotting
            [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
            scene_timerange = scene_timerange + t_start;
            scene_rocket_statematrix = [scene_timerange, scene_statematrix(:,1:6)];
            plot3(scene_rocket_statematrix(:,2),scene_rocket_statematrix(:,3),scene_rocket_statematrix(:,4))
            scatter3(r_target(1),r_target(2),r_target(3),10,"m","filled")
    
            plot3(target_orbit_statematrix(:,2)+state_earth_rendezvous(1),target_orbit_statematrix(:,3)+state_earth_rendezvous(2),target_orbit_statematrix(:,4)+state_earth_rendezvous(3),"k")
            [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(state_earth_rendezvous(1),state_earth_rendezvous(2),state_earth_rendezvous(3),earth_radius,earth_radius,earth_radius,64);
            earth_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[0,0,1]);
    
            station_v_rendezvous = state_earth_rendezvous(4:6) + v_target_offset;
            station_v_rendezvous %what we use to calculate final dv
    
            break %we've arrived
        end
        if toc > 0.25
            break %taking too long
        end
    end
    shooting_good

end

hold on
grid on
axis tight equal
%scatter3(earth_start_statematrix(2),earth_start_statematrix(3),earth_start_statematrix(4),5,"b","filled")
%scatter3(lunar_start_statematrix(2),lunar_start_statematrix(3),lunar_start_statematrix(4),5,"k","filled")
scatter3(rocket_start_statematrix(2),rocket_start_statematrix(3),rocket_start_statematrix(4),10,"r","filled")

lunar_v_norm = lunar_v_initial./norm(lunar_v_initial);
earth_v_norm = earth_v_initial./norm(earth_v_initial);
rocket_start_norm = rocket_start_statematrix(5:7)/norm(rocket_start_statematrix(5:7));

plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b")
plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"k")

% unit_vec = 0.2e8;
% scatter3(target_r_initial(1),target_r_initial(2),target_r_initial(3),5,"m","filled")
% target_v_norm = target_v_initial./norm(target_v_initial);
% quiver3(target_r_initial(1),target_r_initial(2),target_r_initial(3), target_v_norm(1)*unit_vec,target_v_norm(2)*unit_vec,target_v_norm(3)*unit_vec,"m")


%view([45,15])
view([0,90])


function [jacobian, residual_matrix, r_target, state_earth_rendezvous] = fetch_jacobian_and_residual(r_target_offset, ER3BP_timespan, ER3BP_state_inital, lunar_mass, earth_mass, error_tolerance)

    %error_tolerance = 1e-9;

    peturb_factor = 1e-4;

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


function state_out = interp_statematrix_via_TA(statematrix,TA_target)
    row_time = 8;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-TA_target),2);
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], TA_target);
    end
end

function state_out = interp_statematrix_via_time(statematrix,t_target)
    row_time = 1;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-t_target),2);
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], t_target);
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

function state_out = threebody_ODE(t, state_in, earth_mass, lunar_mass)

    G = 6.6743e-11;

    r_rocket = state_in(1:3).';
    v_rocket = state_in(4:6).';
    r_earth = state_in(7:9).';
    v_earth = state_in(10:12).';
    r_luna = state_in(13:15).';
    v_luna = state_in(16:18).';

    r_earthluna = r_earth - r_luna;
    r_unit_earthluna = norm(r_earthluna);

    acc_earth =  -G * lunar_mass * r_earthluna / r_unit_earthluna^3;
    acc_luna = G * earth_mass * r_earthluna / r_unit_earthluna^3;
    
    acc_rocket = (-G*earth_mass)*(r_rocket-r_earth)/(norm(r_rocket-r_earth)^3) + (-G*lunar_mass)*(r_rocket-r_luna)/(norm(r_rocket-r_luna)^3);

    state_out = [v_rocket, acc_rocket, v_earth, acc_earth, v_luna, acc_luna].';
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


function status = stop_too_slow(t,y,flag)
    persistent start_time
    status = 0;
    if strcmp(flag,'init')
        start_time = tic;
    elseif isempty(flag)
        time_max = 0.1; %time limit for cutting off ode solver
        if toc(start_time) > time_max
            status = 1; %stop ode solver
        end
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