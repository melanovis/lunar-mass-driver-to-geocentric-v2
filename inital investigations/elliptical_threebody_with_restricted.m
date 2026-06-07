format compact
clear
clc
clf reset

% ----------

G = 6.6743e-11;
error_tolerance = 1e-12;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6 + 1e6; %m plus atmospheric margin

collision_tolerance = 100e3; %m, the border around objects to assume no collision

lunar_mass = 7.35e22; %kg
earth_mass = 5.972e24;
mu = G*(earth_mass + lunar_mass);

approx_semimajor = 385692.5e3;

pair_eccentricity = 0.0549;

lunar_SOI_radius = approx_semimajor*(1-pair_eccentricity)*(lunar_mass / (3*(earth_mass+lunar_mass)) ) ^ (1/3);

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

earth_init_statematrix = [earth_init_statematrix, earth_TA_list.']; %these are only for finding inital params
lunar_init_statematrix = [lunar_init_statematrix, lunar_TA_list.'];


%% regarding threebody sim
% t_start_norm = 0.075;
% t_start = t_start_norm*earthluna_orbital_period;
% earth_start_statematrix = interp_statematrix_via_time(earth_init_statematrix,t_start);
% lunar_start_statematrix = interp_statematrix_via_time(lunar_init_statematrix,t_start);

TA_start_norm = 0.2;
TA_start = rem(TA_start_norm*360,360);
lunar_start_statematrix = interp_statematrix_via_TA(lunar_init_statematrix,TA_start);
t_start = lunar_start_statematrix(1);
earth_start_statematrix = interp_statematrix_via_time(earth_init_statematrix,t_start);

%rocket_start_statematrix = [t_start, 0, 3.4e8, 0, -1.02e3, 20, 0, 0]; %just test values for now
rocket_start_statematrix = [t_start, lunar_start_statematrix(2:4)+[0,-lunar_radius,0], lunar_start_statematrix(5:7)+[100,-2.35e3,0], 0];

ER3BP_state_inital = [rocket_start_statematrix(2:7), earth_start_statematrix(2:7), lunar_start_statematrix(2:7)];

ER3BP_timespan = [0,earthluna_orbital_period];

[scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
scene_timerange = scene_timerange + t_start;
scene_rocket_statematrix = [scene_timerange, scene_statematrix(:,1:6)];
scene_earth_statematrix = [scene_timerange, scene_statematrix(:,7:12)];
scene_lunar_statematrix = [scene_timerange, scene_statematrix(:,13:18)];

ODE_solve_timeout = false;
if scene_timerange(end) < earthluna_orbital_period + t_start
    ODE_solve_timeout = true;
end

lunar_approach_min = sqrt(sum((scene_rocket_statematrix(:,2:4) - scene_lunar_statematrix(:,2:4)).^2, 2));
%we only want to be checking for collisions once the collision tolerance has been exceeded
ind_lunarcoltol = 1;
while ind_lunarcoltol < numel(lunar_approach_min)
    if lunar_approach_min > collision_tolerance + lunar_radius
        break
    else
        ind_lunarcoltol = ind_lunarcoltol+1;
    end
end
lunar_approach_min_secondary = lunar_approach_min;
lunar_approach_min_secondary(1:ind_lunarcoltol-1) = inf;

[~,ind_min] = min(lunar_approach_min);
lunar_r_closest = norm(scene_rocket_statematrix(ind_min,2:3) - scene_lunar_statematrix(ind_min,2:3)); %used for collision checks
[~,ind_min] = min(lunar_approach_min_secondary);
lunar_r_closest_secondary = norm(scene_rocket_statematrix(ind_min,2:3) - scene_lunar_statematrix(ind_min,2:3)); %used for collision checks


ind_test = ind_min;

earth_approach_min = sqrt(sum((scene_rocket_statematrix(:,2:4) - scene_earth_statematrix(:,2:4)).^2, 2));
[~,ind_min] = min(earth_approach_min);
earth_r_closest = norm(scene_rocket_statematrix(ind_min,2:4) - scene_earth_statematrix(ind_min,2:4));

collision_flag = false;
if lunar_r_closest < lunar_radius || earth_r_closest < earth_radius || lunar_r_closest_secondary < (lunar_radius + collision_tolerance)
    collision_flag = true;
end

%if any high then its not a valid path
check_matrix = [
    collision_flag
    ODE_solve_timeout
];
any(check_matrix)


%subplot(1,2,1)
hold on
grid on
axis equal
scatter3(earth_start_statematrix(2),earth_start_statematrix(3),earth_start_statematrix(4),5,"b","filled")
scatter3(lunar_start_statematrix(2),lunar_start_statematrix(3),lunar_start_statematrix(4),5,"k","filled")
scatter3(rocket_start_statematrix(2),rocket_start_statematrix(3),rocket_start_statematrix(4),10,"r","filled")

scatter3(scene_rocket_statematrix(ind_test,2),scene_rocket_statematrix(ind_test,3),scene_rocket_statematrix(ind_test,4),10,"cyan","filled")
scatter3(scene_lunar_statematrix(ind_test,2),scene_lunar_statematrix(ind_test,3),scene_lunar_statematrix(ind_test,4),10,"cyan","filled")

unit_vec_size = 2e7;
lunar_v_norm = lunar_v_initial./norm(lunar_v_initial);
earth_v_norm = earth_v_initial./norm(earth_v_initial);
rocket_start_norm = rocket_start_statematrix(5:7)/norm(rocket_start_statematrix(5:7));

quiver3(rocket_start_statematrix(2),rocket_start_statematrix(3),rocket_start_statematrix(4),rocket_start_norm(1)*unit_vec_size,rocket_start_norm(2)*unit_vec_size,rocket_start_norm(3)*unit_vec_size,"r")

plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b")
plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"k")
plot3(scene_rocket_statematrix(:,2),scene_rocket_statematrix(:,3),scene_rocket_statematrix(:,4),"r")

view([30,20])


% subplot(1,2,2)
% for n=1:30:height(scene_lunar_statematrix)
% 
%     scatter(nan,nan)
%     hold on
%     grid on
%     axis equal
%     scatter3(scene_rocket_statematrix(n,2),scene_rocket_statematrix(n,3),scene_rocket_statematrix(n,4),10,"b","filled")
%     scatter3(scene_lunar_statematrix(n,2),scene_lunar_statematrix(n,3),scene_lunar_statematrix(n,4),10,"k","filled")
%     plot3(scene_rocket_statematrix(:,2),scene_rocket_statematrix(:,3),scene_rocket_statematrix(:,4),"b")
%     plot3(scene_lunar_statematrix(:,2),scene_lunar_statematrix(:,3),scene_lunar_statematrix(:,4),"k")
% 
%     view([0,90])
%     drawnow()
% 
%     hold off
% end


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

function status = stop_too_slow(t,y,flag)
    persistent start_time
    status = 0;
    if strcmp(flag,'init')
        start_time = tic;
    elseif isempty(flag)
        time_max = 0.2; %time limit for cutting off ode solver
        if toc(start_time) > time_max
            status = 1; %stop ode solver
        end
    end
end