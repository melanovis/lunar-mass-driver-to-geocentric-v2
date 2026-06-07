format compact
clear
clc
clf reset

%----

res = 8;
omega_end_range = linspace(0,360,res+1);
omega_end_range(end) = [];

load("pso_results_collected.mat")

sim_labels = string(fieldnames(pso_results_total));


G = 6.6743e-11;
error_tolerance = 1e-12;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; %m

lunar_obliquity = 6.68;
earth_obliquity = 23.44;

target_altitude = 1371.2e3;
target_radius = target_altitude + earth_radius;

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

timespan = linspace(0,earthluna_orbital_period*5,10e3);

[timerange_init, statematrix_init] = ode113(@(timerange, statematrix) earth_moon_ODE(timerange, statematrix, earth_mass, lunar_mass), timespan, earthlunar_state_intial ,odeset('Reltol',error_tolerance));
earth_init_statematrix = [timerange_init, statematrix_init(:,1:6)];
lunar_init_statematrix = [timerange_init, statematrix_init(:,7:12)];

lunar_peri_state = lunar_init_statematrix(1,:);
earth_apo_state = earth_init_statematrix(1,:);

lunar_TA_list = get_planet_TAlist(lunar_peri_state(2:4),lunar_peri_state(5:7),lunar_init_statematrix);
earth_TA_list = get_planet_TAlist(earth_apo_state(2:4),earth_apo_state(5:7),earth_init_statematrix);
earth_TA_list = rem(earth_TA_list+180,360);

earth_TA_list = accumulate_increasing_TAlist(earth_TA_list);
lunar_TA_list = accumulate_increasing_TAlist(lunar_TA_list);

earth_init_statematrix = [earth_init_statematrix, earth_TA_list.'];
lunar_init_statematrix = [lunar_init_statematrix, lunar_TA_list.'];


for n=1:numel(sim_labels)

    scatter(nan,nan,"k")
    hold on
    grid on
    axis tight equal

    statematrix_spec = [];
    statematrix_spec = pso_results_total.(sim_labels(n)).results_struct.transit_statematrix;

    transit_statematrix = [statematrix_spec(:,1:7), statematrix_spec(:,end)];

    plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b");
    plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"w");

    state_earth_arrival = statematrix_spec(end,8:13);
    state_luna_depart = statematrix_spec(1,14:19);

    [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(state_earth_arrival(1),state_earth_arrival(2),state_earth_arrival(3),earth_radius,earth_radius,earth_radius,64);
    surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0,FaceColor=[0,0,1],FaceAlpha=0.3);
    [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(state_luna_depart(1),state_luna_depart(2),state_luna_depart(3),lunar_radius,lunar_radius,lunar_radius,64);
    surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[repelem(0.5,3)]);

    scatter3(state_luna_depart(1),state_luna_depart(2),state_luna_depart(3),200,"w+");
    text(state_luna_depart(1),state_luna_depart(2),state_luna_depart(3), [""," Luna departure"], Interpreter="latex", FontSize=12,VerticalAlignment="cap",color=[1,1,1]);

    lunar_tick_series = linspace(0,earthluna_orbital_period,36);
    lunar_tick_series(end) = [];
    for m=1:numel(lunar_tick_series)
        [lunar_orbit_tickstate, ~] = interp_statematrix_via_time(lunar_init_statematrix,lunar_tick_series(m));
        lunar_ticks(m) = plot3([lunar_orbit_tickstate(2),lunar_orbit_tickstate(2)],[lunar_orbit_tickstate(3),lunar_orbit_tickstate(3)],[0,lunar_orbit_tickstate(4)],"--",color=[repelem(0.5,3)],LineWidth=0.03);
    end

    patch([transit_statematrix(:,2);nan],[transit_statematrix(:,3);nan],[transit_statematrix(:,4);nan],[transit_statematrix(:,end);nan],'EdgeColor',[1,1,1],linewidth=0.5)

    ind_transitions = find(diff(transit_statematrix(:,end)) ~= 0);
    ind_transitions = [ind_transitions; height(transit_statematrix)];
    burn_list = pso_results_total.(sim_labels(n)).results_struct.burn_dv_series;
    unit_vec_scale = 0.3e8;
    for m=1:height(burn_list)
        desc = string( round(norm(burn_list(m,:))/1e3,3) ) + " km/s burn" + newline + ...
            "t="+round(transit_statematrix(ind_transitions(m),1)./3600,2)+" hours";
        text(transit_statematrix(ind_transitions(m),2),transit_statematrix(ind_transitions(m),3),transit_statematrix(ind_transitions(m),4), [" ", desc], Interpreter="latex", color=[1,0,0],FontSize=10, VerticalAlignment="cap");
        burn_vector = burn_list(m,:);
        burn_vector = burn_vector./norm(burn_vector);
        quiver3( transit_statematrix(ind_transitions(m),2),transit_statematrix(ind_transitions(m),3),transit_statematrix(ind_transitions(m),4), burn_vector(1)*unit_vec_scale,burn_vector(2)*unit_vec_scale,burn_vector(3)*unit_vec_scale, color=[1,0,0],MaxHeadSize=1.2)
        scatter3(transit_statematrix(ind_transitions(m),2),transit_statematrix(ind_transitions(m),3),transit_statematrix(ind_transitions(m),4),20,"r","filled")
        plot3([repelem(transit_statematrix(ind_transitions(m),2),2)],[repelem(transit_statematrix(ind_transitions(m),3),2)],[0,transit_statematrix(ind_transitions(m),4)],"--r",linewidth=0.1)
    end

    target_orbit_stateinital = [pso_results_total.(sim_labels(n)).results_struct.targeting_rv(1,:), pso_results_total.(sim_labels(n)).results_struct.targeting_rv(2,:)] ;
    target_orbit_timespan = [0,sqrt((target_radius^3)*(4*pi^2)/(G*earth_mass))];
    [target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),target_orbit_timespan,target_orbit_stateinital,odeset('Reltol',error_tolerance));
    target_orbit_statematrix = [target_timerange,target_state_matrix];
    for m=2:4
        target_orbit_statematrix(:,m) = + target_orbit_statematrix(:,m) + state_earth_arrival(m-1);
    end
    plot3(target_orbit_statematrix(:,2),target_orbit_statematrix(:,3),target_orbit_statematrix(:,4),"w")

    TA_omega = pso_results_total.(sim_labels(n)).TA_omega;
    launch_params = pso_results_total.(sim_labels(n)).results_struct.MD_params;
    dv_total = pso_results_total.(sim_labels(n)).results_struct.dv_total;

    legend( "MD allocation phase transits - MDR " + erase(sim_labels(n),"ind_s_") + newline + ...
        "$\nu_l$ = "+ TA_omega(1)+"$^{\circ}$, "+"$\Omega_t$ = "+ TA_omega(2)+"$^{\circ}$" + newline + ...
        "launch v = " + round(launch_params(1)/1e3,3) + " km/s" + ", rocket $\Delta$v = " + round(dv_total/1e3,3)+" km/s" + newline + ...
        "launch $\phi$ = " + round(launch_params(3),2)+"$^{\circ}$" + ", launch $\lambda$ = " + round(launch_params(2),2)+"$^{\circ}$" + newline + ...
        "launch $\theta$ = " +round(launch_params(4),2)+"$^{\circ}$"+ ", launch $\epsilon$ = "+round(launch_params(5),6)+"$^{\circ}$", ...
        textcolor="w", FontSize=18, Interpreter="latex",location="northeast")
    legend boxoff

    view([30,10])
    
    ax = gca;
    ax.FontSize = 20;

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
    
    zlim("padded")

    set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

    hold off
    drawnow()

    im_raw = getframe(gcf);
    imwrite(im_raw.cdata, "MD_allocation_"+sim_labels(n)+".png");
    fprintf("-\n")
end


function [state_out, inds_between] = interp_statematrix_via_time(statematrix,t_target)
    row_time = 1;
    [~,inds_between] = mink(abs(statematrix(:,row_time)-t_target),2);
    inds_between = sort(inds_between);
    for n=1:width(statematrix)
        state_out(n) = interp1( [statematrix(inds_between(1),row_time),statematrix(inds_between(2),row_time)], [statematrix(inds_between(1),n),statematrix(inds_between(2),n)], t_target,"linear","extrap");
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