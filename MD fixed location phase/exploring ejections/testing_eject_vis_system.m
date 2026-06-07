format compact
clear
clc
clf reset

% ----------

load("craters_formatted.mat")

eject_v = 2.33e3;

TA_start_norm = 0.99;

launch_long = 0;
launch_lat = 0;
% launch_long = -45;
% launch_lat = 45;
launch_azimuth = 90;
launch_elevation = 0;

launch_azimuth = rem(launch_azimuth,360);

earth_SOI_approx = 1.5e9;

G = 6.6743e-11;
error_tolerance = 1e-13;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; 

collision_tolerance = 100e3; %m, the border around objects to assume no collision

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
    launch_lat = sign(launch_lat)*(90-1e-3);
end


launch_r = (cosd(launch_lat)*cosd(launch_long))*null_norm + (cosd(launch_lat)*sind(launch_long))*equator_90_norm + sind(launch_lat)*pole_norm;

launch_r_norm = (launch_r./norm(launch_r));
launch_r_lunar = launch_r_norm * (lunar_radius+10);

launch_tang_plane = pole_norm - dot(pole_norm,launch_r_norm)*launch_r_norm;
launch_ref_north = launch_tang_plane./norm(launch_tang_plane);
launch_ref_east = cross(launch_r_norm, launch_ref_north);

tang_vector = cosd(-launch_azimuth)*launch_ref_north + sind(-launch_azimuth)*launch_ref_east;
eject_vector = cosd(launch_elevation)*tang_vector + sind(launch_elevation)*launch_r_norm;

eject_vector_norm = eject_vector/norm(eject_vector);

eject_v_geo = eject_vector_norm*eject_v + lunar_start_statematrix(5:7);
eject_r_geo = lunar_start_statematrix(2:4) + launch_r_lunar;

total_transit_statematrix = [];

scene_time = earthluna_orbital_period*2.5;

scene_timerange = [];
scene_statematrix = [];
rocket_start_statematrix = [t_start, eject_r_geo, eject_v_geo];
ER3BP_state_inital = [rocket_start_statematrix(2:7), earth_start_statematrix(2:7), lunar_start_statematrix(2:7)];
ER3BP_timespan = [0,scene_time];
[scene_timerange, scene_statematrix] = ode113(@(scene_timerange, scene_statematrix) threebody_ODE(scene_timerange, scene_statematrix, earth_mass, lunar_mass), ER3BP_timespan, ER3BP_state_inital ,odeset('Reltol',error_tolerance, 'OutputFcn', @stop_too_slow));
scene_timerange = scene_timerange + t_start;
scene_eject_statematrix = [scene_timerange, scene_statematrix]; %statematrix before first burn

total_transit_statematrix = scene_eject_statematrix;


lunar_approach = sqrt(sum((total_transit_statematrix(:,2:4) - total_transit_statematrix(:,14:16)).^2, 2));
%we only want to be checking for collisions once the collision tolerance has been exceeded
ind_lunarcoltol = 1;
while ind_lunarcoltol < numel(lunar_approach)
    if lunar_approach(ind_lunarcoltol) > collision_tolerance + lunar_radius
        break
    else
        ind_lunarcoltol = ind_lunarcoltol+1;
    end
end
lunar_approach_secondary = lunar_approach;
lunar_groundingcheck = lunar_approach_secondary(1:ind_lunarcoltol-1);
lunar_approach_secondary(1:ind_lunarcoltol-1) = inf;

[~,ind_min_l1] = min(lunar_approach);
lunar_r_closest = norm(total_transit_statematrix(ind_min_l1,2:4) - total_transit_statematrix(ind_min_l1,14:16));

[~,ind_min_l2] = min(lunar_approach_secondary);
lunar_r_closest_secondary = norm(total_transit_statematrix(ind_min_l2,2:4) - total_transit_statematrix(ind_min_l2,14:16)); 

earth_approach = sqrt(sum((total_transit_statematrix(:,2:4) - total_transit_statematrix(:,8:10)).^2, 2));
[~,ind_min_e] = min(earth_approach);
earth_r_closest = norm(total_transit_statematrix(ind_min_e,2:4) - total_transit_statematrix(ind_min_e,8:10));

collision_flag = false;
if lunar_r_closest < lunar_radius || earth_r_closest < (earth_radius + collision_tolerance) || lunar_r_closest_secondary < (lunar_radius + collision_tolerance) || any(lunar_groundingcheck < lunar_radius)
    collision_flag = true;
end

worst_collision_intersect = 0;
if collision_flag
    worst_collision_intersect = min([lunar_r_closest - lunar_radius, lunar_r_closest_secondary - (lunar_radius + collision_tolerance), earth_r_closest - (earth_radius + collision_tolerance) ]);
end

transit_time = total_transit_statematrix(end,1) - total_transit_statematrix(1,1);

check_matrix = [
~collision_flag
~any(any(isnan(total_transit_statematrix)))
];

height_original = height(total_transit_statematrix);

if collision_flag
    if lunar_r_closest < lunar_radius
        ind_collision = ind_min_l1;
        collision_str = "lunar collision!";
    elseif lunar_r_closest_secondary < (lunar_radius + collision_tolerance)
        ind_collision = ind_min_l2;
        collision_str = "lunar collision!";
    else 
        ind_collision = ind_min_e;
        collision_str = "earth collision!";
    end

    total_transit_statematrix(ind_collision+1:end,:)=[];
end

deep_space_flag = false;
if max(earth_approach) > earth_SOI_approx && ~collision_flag
    deep_space_flag = true;
    deep_str = "carried away onto helio orbit!";
    [~,ind_SOI_exit] = min(abs(earth_approach-earth_SOI_approx));
    total_transit_statematrix(ind_SOI_exit+1:end,:) = [];
end


cmap = interp1([0,0.01,0.2,0.4,0.6,0.8,1], [[repelem(0,3)]; [repelem(0.1,3)]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));

clf reset

scatter(nan,nan,"k");
%set(gcf, 'Renderer', 'painters')

for n=1:2
    
    if n==1
        subplot(2,2,1)
    else
        subplot(2,2,[2,4])
        scatter(nan,nan,"k");
    end

    hold on
    grid on
    axis tight equal
    scatter3(lunar_start_statematrix(2),lunar_start_statematrix(3),lunar_start_statematrix(4),5,"w","filled")
    plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b")
    plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"w")
    
    ax = gca;
    ax.FontSize = 20;

    ind_test = 1;
    state_earth_arrival = total_transit_statematrix(end,8:13);
    state_luna_depart = total_transit_statematrix(1,14:19);
    state_rocket_test = total_transit_statematrix(1,2:4);
    
    [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(state_luna_depart(1),state_luna_depart(2),state_luna_depart(3),lunar_radius,lunar_radius,lunar_radius,64);
    luna_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[repelem(0.7,3)]);
    luna_marker = scatter3(state_luna_depart(1),state_luna_depart(2),state_luna_depart(3),300,"wx");
    
    scatter3(state_rocket_test(1),state_rocket_test(2),state_rocket_test(3),5,"w","filled")
    
    %plot3(total_transit_statematrix(:,2),total_transit_statematrix(:,3),total_transit_statematrix(:,4),"-r")
    % trail_patch = patch([total_transit_statematrix(:,2); nan],[total_transit_statematrix(:,3); nan],[total_transit_statematrix(:,4); nan],'w','EdgeColor','w',...
    %     'FaceVertexAlphaData',linspace(1, (height_original - height(total_transit_statematrix)) / height_original , height(total_transit_statematrix)+1).','AlphaDataMapping','none',...
    %     'EdgeAlpha','interp', 'LineWidth', 2);
    trail_patch = patch([total_transit_statematrix(:,2); nan],[total_transit_statematrix(:,3); nan],[total_transit_statematrix(:,4); nan], [linspace(0,abs(total_transit_statematrix(1,1)-total_transit_statematrix(end,1))/(3600*24),height(total_transit_statematrix)).'; nan],'EdgeColor','interp', 'LineWidth', 2);
    colormap(flip(cmap));
    clim([0,scene_time/(3600*24)])

    if n==2
        h = colorbar;
        set(get(h,'label'),'string','time (days)', Interpreter="latex", FontSize=20);
        h.Color = [1,1,1];
        h.Location = "southoutside";
    end

    if collision_flag
        scatter3(total_transit_statematrix(end,2),total_transit_statematrix(end,3),total_transit_statematrix(end,4),10,"r","filled")
        text(total_transit_statematrix(end,2),total_transit_statematrix(end,3),total_transit_statematrix(end,4), "$\enspace$" + collision_str, Interpreter="latex", FontSize=14, Color = [1,0,0])
    end

    if n==1
        view([-90,0])
    else
        view([0,90])
    end
    zlim("padded")

    if n==2
        xlim([-1,1].*6e8)
        ylim([-1,1].*6e8)
        zlim([-1,1].*20e8)
    else
        xlim([-1,1].*20e8)
        ylim([-1,1].*5.5e8)
        zlim([-1,1].*3.1e8)
    end

    set(gcf, 'Color', [0,0,0])
    set(ax, 'Color', [0,0,0])
    
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.ZColor = 'w';
    ax.Title.Color  = 'w';
    ax.XLabel.Color = 'w';
    ax.YLabel.Color = 'w';
    
    ax.GridColor = [1 1 1];
    ax.MinorGridColor = [1 1 1];
    ax.GridAlpha = 0.2;

    if n==2
        legend_str = "eject v = " + round(eject_v/1e3,4) + " km/s" + newline +"eject lunar $\nu$ = "+round(rem(TA_start_norm*360,360),2) + "$^{\circ}$" + newline + ...
        "launch $\theta$ = " + round(launch_azimuth,3) +"$^{\circ}$" + ", " + "launch $\epsilon$ = " + round(launch_elevation,3) +"$^{\circ}$" + newline +...
        "launch $\phi$ = " + round(launch_lat,2)+"$^{\circ}$" + ", launch $\lambda$ = " + round(launch_long,2)+"$^{\circ}$";

        legend(legend_str,textcolor="w", FontSize=15, Interpreter="latex",location="northwest")
        legend boxoff
    end

end

subplot(2,2,3)
hold on
grid on
axis tight equal
xlim([-180,180])
ylim([-90,90])
surfacemap = imread('lunar_map.png');
surfacemap = flipud(surfacemap);

wm = image(surfacemap,'xdata',[-180,180],'ydata',[-90 90]);

for n=1:height(crater_labels)
    if crater_latlong(n,2) < 150 && crater_diameter(n) > 210
        scatter(crater_latlong(n,2),crater_latlong(n,1),"rx")
        text(crater_latlong(n,2),crater_latlong(n,1)," "+string(crater_labels(n)),fontsize = 10,Color=[1,1,1])
    end
end

eject_vector = [
cosd(launch_azimuth), sind(launch_azimuth)
-sind(launch_azimuth), cosd(launch_azimuth)
]*[0;1];

unit_vec_size = interp1([0,3e3],[0,40],eject_v,"linear","extrap");
scatter(launch_long,launch_lat,150,"w+")
scatter(launch_long,launch_lat,30,"r","filled")
quiver(launch_long,launch_lat,eject_vector(1)*unit_vec_size,eject_vector(2)*unit_vec_size,"r",MaxHeadSize=1)

xlabel("longitude", Interpreter="latex", FontSize=20)
ylabel("latitude", Interpreter="latex", FontSize=20)

ax = gca;
ax.FontSize = 20;
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';
ax.Title.Color  = 'w';
ax.XLabel.Color = 'w';
ax.YLabel.Color = 'w';

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

hold off
drawnow()




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
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        if statematrix(inds_between(1),n) == statematrix(inds_between(2),n)
            unique_term = 1e-12;
        else
            unique_term = 0;
        end
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time) + unique_term], [statematrix(inds_between(1),n),statematrix(inds_between(2),n) + unique_term], t_target);
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

