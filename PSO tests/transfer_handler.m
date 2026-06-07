function [return_struct,target_series] = transfer_handler(input_vector, TA_start_norm, target_radius, target_omega, earth_init_statematrix, lunar_init_statematrix, earthluna_orbital_period)

warning('off', 'all');

max_manuevers = 2;

manuevers = round(interp1([0,1],[1,max_manuevers], input_vector(1) ));

eject_v = interp1([0,1],[1e3,3e3], input_vector(2) ); %m/s

dt_final = input_vector(3); %normalised
dt_final_shootinguess = input_vector(4); %normalised

target_TA = interp1([0,1],[0,360], input_vector(5) );

%launch params, will be fixed later
launch_long = interp1([0,1],[-90,90], input_vector(6) ); %lambda
launch_lat = interp1([0,1],[-90,90], input_vector(7) ); %phi
launch_azimuth = interp1([0,1],[0,360], input_vector(8) ); %theta, always defined relative to north heading, moves around clockwise
launch_elevation = interp1([0,1],[0,30], input_vector(9) ); %epsilon, relative to surface tangent

earth_SOI_approx = 1.5e9;

manuever_matrix = input_vector(10:10 + max_manuevers*6 -1 );
manuever_matrix = reshape(manuever_matrix,max_manuevers,[]);

G = 6.6743e-11;
error_tolerance = 1e-10;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; 

collision_tolerance = 150e3; %m, the border around objects to assume no collision

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


TA_start = rem(TA_start_norm*360,360);
lunar_start_statematrix = interp_statematrix_via_TA(lunar_init_statematrix,TA_start);
t_start = lunar_start_statematrix(1);
earth_start_statematrix = interp_statematrix_via_time(earth_init_statematrix,t_start);

%% regarding launch

perp_lunar_plane = cross(lunar_start_statematrix(2:4),[0,0,1]);
perp_lunar_plane = perp_lunar_plane./norm(perp_lunar_plane);

pole_unshifted = [0,0,1];
pole_theta = lunar_inc - lunar_obliquity;
pole_corrected = pole_unshifted*cosd(pole_theta) + cross(perp_lunar_plane,pole_unshifted)*sind(pole_theta);

null_point = pole_corrected*cosd(90) + cross(perp_lunar_plane,pole_corrected)*sind(90); %this is 0,0 long/lat
equator_90 = cross(pole_corrected,null_point);

pole_norm = pole_corrected./norm(pole_corrected);
null_norm = null_point./norm(null_point);
equator_90_norm = equator_90./norm(equator_90);

if abs(launch_lat) == 90
    launch_lat = sign(launch_lat)*(90-1e-4);
end

launch_r = (cosd(launch_lat)*cosd(launch_long))*null_norm + (cosd(launch_lat)*sind(launch_long))*equator_90_norm + sind(launch_lat)*pole_norm;
launch_r_norm = (launch_r./norm(launch_r));
launch_r_lunar = launch_r_norm * (lunar_radius+1);

launch_tang_plane = pole_norm - dot(pole_norm,launch_r_norm)*launch_r_norm;
launch_ref_north = launch_tang_plane./norm(launch_tang_plane);
launch_ref_east = cross(launch_r_norm, launch_ref_north);

tang_vector = cosd(-launch_azimuth)*launch_ref_north + sind(-launch_azimuth)*launch_ref_east;
eject_vector = cosd(launch_elevation)*tang_vector + sind(launch_elevation)*launch_r_norm;
eject_vector_norm = eject_vector./norm(eject_vector);

eject_v_geo = eject_vector_norm*eject_v + lunar_start_statematrix(5:7);
eject_r_geo = lunar_start_statematrix(2:4) + launch_r_lunar;


%% building target orbit

target_inc = earth_obliquity + lunar_inc;

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

%% manuever sequencing

target_series = [];
dt_lambert_series = [];
dt_burn_series = [];
throttle_series = [];
for n=1:manuevers

    target_r = interp1([0,1],[earth_radius+800e3, earth_SOI_approx/2],manuever_matrix(n,2));
    target_theta = interp1([0,1],[0,360],manuever_matrix(n,3));
    target_h = interp1([0.1,1],[-1e8, 1e8],manuever_matrix(n,4));

    dt_burn_series(n,1) = interp1([0,1],[0,1-1e-3],manuever_matrix(n,1));
    dt_lambert_series(n,1) = manuever_matrix(n,5);

    if manuever_matrix(n,6) > 0.9
        throttle_series(n,1) = 1;
    else
        throttle_series(n,1) = 10 ^ interp1([0,0.9],[-5, 0.5],manuever_matrix(n,6));
    end
    
    r_point = [target_r*cosd(target_theta), target_r*sind(target_theta), target_h];
    target_series(n,:) = r_point;
end
% target_series
% dt_lambert_series


scene_timerange = [];
scene_statematrix = [];
rocket_start_statematrix = [t_start, eject_r_geo, eject_v_geo];
ER3BP_state_inital = [rocket_start_statematrix(2:7), earth_start_statematrix(2:7), lunar_start_statematrix(2:7)];
ER3BP_timespan = [0,earthluna_orbital_period];
[scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
scene_timerange = scene_timerange + t_start;
scene_eject_statematrix = [scene_timerange, scene_statematrix]; %statematrix before first burn

total_transit_statematrix = [];

burn_r_series = [];
burn_dv_series = []; 
ODE_solve_stable = true;

if manuevers ~= 0

    statematrix_prev = scene_eject_statematrix;
    
    leg_statematrix = [];
        
    for ind_m = 1:manuevers
       
        burn_time = interp1([0,1],[ min(statematrix_prev(:,1)), max(statematrix_prev(:,1)) ],dt_burn_series(ind_m,1));
        [state_at_burn, inds_between] = interp_statematrix_via_time(statematrix_prev,burn_time);
        leg_statematrix = [statematrix_prev(1:min(inds_between),:); state_at_burn];
        leg_statematrix = [leg_statematrix, repelem(ind_m,height(leg_statematrix)).'];
    
        total_transit_statematrix = [total_transit_statematrix; leg_statematrix]; %add on index of leg
    
        [v_1,~,~] = lambert_solve_anticlockwise(state_at_burn(2:4), [target_series(ind_m,:)], dt_lambert_series(ind_m,1), G*earth_mass);
    
        t_start = state_at_burn(1);
        new_path_start = [t_start, state_at_burn(2:4), v_1.*throttle_series(ind_m,1)];
    
        ER3BP_state_inital = [new_path_start(2:7), state_at_burn(8:13), state_at_burn(14:19)]; %new rocket, earth and luna states
        if any(isnan(new_path_start))
            ER3BP_state_inital = leg_statematrix(end-1,2:end-1);
            t_start = leg_statematrix(end-1,1);
            ODE_solve_stable = false;
        end
        ER3BP_timespan = [0,earthluna_orbital_period];
        [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
        scene_timerange = scene_timerange + t_start;
        new_path_statematrix = [scene_timerange, scene_statematrix];
    
        burn_vector = new_path_statematrix(1,5:7) - state_at_burn(5:7);
        burn_r_series = [burn_r_series; state_at_burn(2:4)];
        burn_dv_series = [burn_dv_series; burn_vector];
    
        statematrix_prev = new_path_statematrix;
    end
else
    new_path_statematrix = scene_eject_statematrix;
    ind_m = 0;
    if any(any(isnan(new_path_statematrix)))
        ODE_solve_stable = false;
    end
end

burn_time = interp1([0,1],[ min(new_path_statematrix(:,1)), max(new_path_statematrix(:,1)) ],dt_final);
[state_at_burn, inds_between] = interp_statematrix_via_time(new_path_statematrix,burn_time);
penultimate_leg_statematrix = [new_path_statematrix(1:min(inds_between),:); state_at_burn];
penultimate_leg_statematrix = [penultimate_leg_statematrix, repelem(ind_m+1,height(penultimate_leg_statematrix)).'];

total_transit_statematrix = [total_transit_statematrix; penultimate_leg_statematrix];

%% shooting method implementation

t_startshoot = total_transit_statematrix(end,1);
r_shooting_start = total_transit_statematrix(end,2:4);
v_shooting_start = total_transit_statematrix(end,5:7);
earth_startshoot_statematrix = total_transit_statematrix(end,8:13);
lunar_startshoot_statematrix = total_transit_statematrix(end,14:19);

if isnan(t_startshoot)
    t_startshoot = total_transit_statematrix(end-1,1);
    r_shooting_start = total_transit_statematrix(end-1,2:4);
    v_shooting_start = total_transit_statematrix(end-1,5:7);
    earth_startshoot_statematrix = total_transit_statematrix(end-1,8:13);
    lunar_startshoot_statematrix = total_transit_statematrix(end-1,14:19);
    ODE_solve_stable = false;
end

t_startshoot = double(t_startshoot);
r_shooting_start = double(r_shooting_start);
v_shooting_start = double(v_shooting_start);
earth_startshoot_statematrix = double(earth_startshoot_statematrix);
lunar_startshoot_statematrix = double(lunar_startshoot_statematrix);

[v_1, ~, delta_t] = lambert_solve_anticlockwise(r_shooting_start, [0, -target_radius ,0], dt_final_shootinguess, G*earth_mass);
v_guess_inital = v_1;

rocket_startshoot_statematrix = [t_startshoot, r_shooting_start, v_guess_inital];
ER3BP_state_inital = [rocket_startshoot_statematrix(2:7), earth_startshoot_statematrix, lunar_startshoot_statematrix];
[scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
scene_timerange = scene_timerange + t_startshoot;
scene_rocket_statematrix = [scene_timerange, scene_statematrix(:,1:6)];

ER3BP_timespan = double(t_startshoot + [0, delta_t]);

shooting_good = false;
lastleg_statematrix = [];

tic
for n=1:30

    [jacobian, residual_matrix, r_target, state_earth_rendezvous] = fetch_jacobian_and_residual(r_target_offset, ER3BP_timespan, ER3BP_state_inital, lunar_mass, earth_mass, error_tolerance);
    
    [u_mat,s_mat,v_mat] = svd(jacobian);
    step_matrix = -v_mat * (s_mat \ (u_mat' * residual_matrix));
    v_step = step_matrix(1:3);
    t_step = step_matrix(4);

    ER3BP_state_inital(5:7) = ER3BP_state_inital(5:7) + v_step.';
    ER3BP_timespan(2) = ER3BP_timespan(2) + t_step;

    target_r_rendezvous = state_earth_rendezvous(1:3) + r_target_offset;

    if norm(residual_matrix) < 10
        shooting_good = true;
        
        [scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
        lastleg_statematrix = [scene_timerange, scene_statematrix];
        lastleg_statematrix = [lastleg_statematrix, repelem(ind_m+2,height(lastleg_statematrix)).'];

        parking_v_rendezvous = state_earth_rendezvous(4:6) + v_target_offset;

        burn_vector_lastleg_start = lastleg_statematrix(1,5:7) - v_shooting_start;
        burn_vector_lastleg_end = parking_v_rendezvous - lastleg_statematrix(end,5:7);
        burn_r_series = [burn_r_series; lastleg_statematrix(1,2:4); lastleg_statematrix(end,2:4)];
        burn_dv_series = [burn_dv_series; burn_vector_lastleg_start; burn_vector_lastleg_end];

        break %we've arrived
    end
    if toc > 0.325
        break %taking too long
    end
end

total_transit_statematrix = [total_transit_statematrix; lastleg_statematrix];

%checking transitions
ind_transitions = find(diff(total_transit_statematrix(:,end)) ~= 0);
for n=1:numel(ind_transitions)
    transition_delta = total_transit_statematrix(ind_transitions(n)+1,1:4) - total_transit_statematrix(ind_transitions(n),1:4);
    transition_check(n) = true;
    %if transition_delta 
    if all(~isnan(transition_delta))
        if numel(unique(logical(transition_delta))) ~= 1
            transition_check(n) = false;
        end
    end
end

lunar_approach = sqrt(sum((total_transit_statematrix(:,2:4) - total_transit_statematrix(:,14:16)).^2, 2));
%we only want to be checking for collisions once the collision tolerance has been exceeded
ind_lunarcoltol = 1;
while ind_lunarcoltol < numel(lunar_approach)
    if lunar_approach > collision_tolerance + lunar_radius
        break
    else
        ind_lunarcoltol = ind_lunarcoltol+1;
    end
end
lunar_approach_secondary = lunar_approach;
lunar_approach_secondary(1:ind_lunarcoltol-1) = inf;

[~,ind_min] = min(lunar_approach);
lunar_r_closest = norm(total_transit_statematrix(ind_min,2:4) - total_transit_statematrix(ind_min,14:16));
[~,ind_min] = min(lunar_approach_secondary);
lunar_r_closest_secondary = norm(total_transit_statematrix(ind_min,2:4) - total_transit_statematrix(ind_min,14:16)); 

earth_approach = sqrt(sum((total_transit_statematrix(:,2:4) - total_transit_statematrix(:,8:10)).^2, 2));
[~,ind_min] = min(earth_approach);
earth_r_closest = norm(total_transit_statematrix(ind_min,2:4) - total_transit_statematrix(ind_min,8:10));

collision_flag = false;
if lunar_r_closest < lunar_radius || earth_r_closest < (earth_radius + collision_tolerance) || lunar_r_closest_secondary < (lunar_radius + collision_tolerance)
    collision_flag = true;
end

worst_collision_intersect = 0;
if collision_flag
    worst_collision_intersect = min([lunar_r_closest - lunar_radius, lunar_r_closest_secondary - (lunar_radius + collision_tolerance), earth_r_closest < (earth_radius + collision_tolerance)]);
end


%% thrust calculation system
impulsive_assumption_angle = 3;
acc_req_series = [];
burn_indexes = find(diff(total_transit_statematrix(:,end)) ~= 0);
thrust_check_good = true;
for n=1:height(burn_indexes)
    ind_centroid = burn_indexes(n);
    %fprintf("---\n")
    r_burn_centroid = total_transit_statematrix(ind_centroid,2:4);

    theta = 0;
    ind_spec = ind_centroid;
    while theta < impulsive_assumption_angle/2
        r_spec = total_transit_statematrix(ind_spec,2:4) - r_burn_centroid;
        theta = acosd( dot(r_burn_centroid,r_spec) / (norm(r_burn_centroid)*norm(r_spec)) );
        ind_spec = ind_spec+1;
        if ind_spec > height(total_transit_statematrix)
            thrust_check_good = false;
            ind_spec = height(total_transit_statematrix);
            break
        end
    end
    thrust_time_forward = total_transit_statematrix(ind_spec,1);

    theta = 0;
    ind_spec = ind_centroid;
    while theta < impulsive_assumption_angle/2
        r_spec = total_transit_statematrix(ind_spec,2:4) - r_burn_centroid;
        theta = acosd( dot(r_burn_centroid,r_spec) / (norm(r_burn_centroid)*norm(r_spec)) );
        ind_spec = ind_spec-1;
        if ind_spec < 1
            thrust_check_good = false;
            ind_spec = 1;
            break
        end
    end
    thrust_time_backward = total_transit_statematrix(ind_spec,1);

    thrust_period_time = abs(thrust_time_forward - thrust_time_backward);

    burn_acc_req = norm(burn_dv_series(n,:))/thrust_period_time;
    acc_req_series = [acc_req_series; burn_acc_req];
end
%accounting for final burn
theta = 0;
r_transit_end = total_transit_statematrix(end,2:4);
ind_spec = height(total_transit_statematrix);
while theta < impulsive_assumption_angle
    r_spec = total_transit_statematrix(ind_spec,2:4) - r_transit_end;
    theta = acosd( dot(r_transit_end,r_spec) / (norm(r_transit_end)*norm(r_spec)) );
    ind_spec = ind_spec-1;
    if ind_spec < 1
        ind_spec = 1;
        thrust_check_good = false;
        break
    end
end

if ~isempty(burn_dv_series)
    thrust_period_time = abs(total_transit_statematrix(end,1) - total_transit_statematrix(ind_spec,1));
    burn_acc_req = norm(burn_dv_series(end,:))/thrust_period_time;
    acc_req_series = [acc_req_series; burn_acc_req];
else
    thrust_check_good = false;
end

if any(isnan(acc_req_series))
    thrust_check_good = false;
end
if ~shooting_good && manuevers == 0
    acc_req_series = [nan];
    burn_dv_series = [nan,nan,nan]; 
    burn_r_series = [nan,nan,nan];
end


dv_total = 0;
for n=1:height(burn_dv_series)
    dv_spec = norm( burn_dv_series(n,:) );
    if isnan(dv_spec)
        dv_spec = 20e3;
    end
    dv_total = dv_total + dv_spec;
end
%dv_total


collision_penalty = exp( -max([log10(abs(worst_collision_intersect)+1e-11), 0])/20 ); %fitness component for collision

deep_space_penalty = 1;
if max(lunar_approach) > earth_SOI_approx/1.5
    %out too far!
    deep_space_penalty = exp( -log10(abs(max(lunar_approach) - earth_SOI_approx/1.5)) / 20 );
end

transit_time = total_transit_statematrix(end,1) - total_transit_statematrix(1,1);

if shooting_good
    thrust_penalty = exp( -max(acc_req_series)/1e3 ); %fitness component for thrust 
    time_penality = exp( -log10(transit_time)/1e3 ); %fitness component for total transit time
    shooting_penalty = 1.5;
else
    thrust_penalty = 0.8;
    time_penality = 0.8;
    shooting_penalty = exp( -log10( norm( total_transit_statematrix(end,2:4) - target_r_rendezvous ) )/70 ); %if shooting bad, we can bring the end closer to the target
end

if isnan(thrust_penalty)
    thrust_penalty = 0.8;
end
if isnan(shooting_penalty)
    shooting_penalty = 0.1;
end

check_matrix = [
ODE_solve_stable
shooting_good
~collision_flag
thrust_check_good
all(~isnan(acc_req_series))
~any(any(isnan(total_transit_statematrix)))
all(transition_check)
];


absent_lastleg_modifier = 0;
if any(~check_matrix)
    absent_lastleg_modifier = sum(~check_matrix)*2e4;
end

eject_velocity_penalty = exp( -eject_v/1e5);

fitness = 5e5/(dv_total + (absent_lastleg_modifier) );

fitness = fitness * shooting_penalty * time_penality * collision_penalty * eject_velocity_penalty * deep_space_penalty * thrust_penalty;

% if transit_time > 92*3600*24 %92 day limit
%     fitness = fitness/10;
% end

if max(acc_req_series) < 6*9.81
    fitness = fitness*1.5;
end

%fitness

return_struct = struct();
return_struct = setfield(return_struct,"fitness",fitness);
return_struct = setfield(return_struct,"acc_series",acc_req_series);
return_struct = setfield(return_struct,"transit_statematrix",total_transit_statematrix);
return_struct = setfield(return_struct,"check_matrix",check_matrix);
return_struct = setfield(return_struct,"burn_dv_series",burn_dv_series);
return_struct = setfield(return_struct,"burn_r_series",burn_r_series);
return_struct = setfield(return_struct,"dv_total",dv_total);
return_struct = setfield(return_struct,"dt_total",transit_time);

MD_param_return = [
eject_v
launch_long
launch_lat
launch_azimuth
launch_elevation
];

return_struct = setfield(return_struct,"MD_params",MD_param_return);

targeting_rv = [
r_target_offset
v_target_offset
];

return_struct = setfield(return_struct,"targeting_rv",targeting_rv);

if all(check_matrix)
    fprintf("#")
else
    fprintf("-")
end



% cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3 ));
% scatter(nan,nan)
% hold on
% grid on
% axis tight equal
% scatter3(lunar_start_statematrix(2),lunar_start_statematrix(3),lunar_start_statematrix(4),5,"k","filled")
% plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b")
% plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"k")
% 
% [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(lunar_start_statematrix(2),lunar_start_statematrix(3),lunar_start_statematrix(4),lunar_radius,lunar_radius,lunar_radius,64);
% lunar_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[repelem(0.5,3)]);
% 
% patch([total_transit_statematrix(:,2);nan],[total_transit_statematrix(:,3);nan],[total_transit_statematrix(:,4);nan],[total_transit_statematrix(:,end);nan],'EdgeColor','interp',linewidth=1.5)
% colormap(cmap)
% clim([1,max(total_transit_statematrix(:,end))+1])
% 
% scatter3(r_shooting_start(1),r_shooting_start(2),r_shooting_start(3),10,"m","filled")
% 
% unitvec = 0.5e8;
% for n=1:height(burn_r_series)
%     burn_norm = burn_dv_series(n,:)./norm(burn_dv_series(n,:));
%     scatter3(burn_r_series(n,1),burn_r_series(n,2),burn_r_series(n,3),10,"r","filled")
%     text(burn_r_series(n,1),burn_r_series(n,2),burn_r_series(n,3), string(norm(burn_dv_series(n,:))./1e3)+" km/s manuever" )
%     quiver3(burn_r_series(n,1),burn_r_series(n,2),burn_r_series(n,3),burn_norm(1)*unitvec,burn_norm(2)*unitvec,burn_norm(3)*unitvec,"r")
% end
% 
% view([30,30])
% 
% set(gcf, 'Color', [1,1,1])
% set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
% drawnow()
% hold off

% xlim(lunar_start_statematrix(2)+[-1e7,1e7])
% ylim(lunar_start_statematrix(3)+[-1e7,1e7])
% zlim(lunar_start_statematrix(4)+[-1e7,1e7])



end



function state_out = interp_statematrix_via_TA(statematrix,TA_target)
    row_time = 8;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-TA_target),2);
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], TA_target);
    end
end

function [state_out, inds_between] = interp_statematrix_via_time(statematrix,t_target)
    row_time = 1;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-t_target),2);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], t_target);
    end
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
