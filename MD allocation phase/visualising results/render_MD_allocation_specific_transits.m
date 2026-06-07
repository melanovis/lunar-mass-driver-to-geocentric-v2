format compact
clear
clc
clf reset

%----

fps = 60;
frames = 40 * 60;

record_video = true;

yaw_bounds = [0,360]+45;
pitch_bounds = [12,12];
yaw_revolutions = 0.25;
pitch_revolutions = 1;

load("pso_results_collected.mat")

sim_labels = string(fieldnames(pso_results_total));

res = 8;
omega_end_range = linspace(0,360,res+1);
omega_end_range(end) = [];

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


for omega_index = 1:numel(omega_end_range)

    clf reset

    %omega_index = 1;
    
    view_type = 2;
    
    %finding boundaries
    time_bounds = [];
    bounds = zeros(7,2);
    time_ending_series = [];
    for n=1:numel(sim_labels)
        statematrix_spec = pso_results_total.(sim_labels(n)).results_struct.transit_statematrix;
        if ismember(pso_results_total.(sim_labels(n)).TA_omega(2), omega_end_range(omega_index))
            for m = 1:4
                bounds(m,:) = [min(min(statematrix_spec(:,m)), bounds(m,1)), max( max(statematrix_spec(:,m)), bounds(m,2) )];
            end
            time_ending_series = [time_ending_series, statematrix_spec(end,1)];
        end
    end
    for m = 2:4
        bounds(m,:) = [min(min(lunar_init_statematrix(:,m)), bounds(m,1)), max( max(lunar_init_statematrix(:,m)), bounds(m,2) )];
    end
    x_bounds = bounds(2,:);
    y_bounds = bounds(3,:);
    z_bounds = bounds(4,:);
    time_bounds = bounds(1,:);
    time_bounds(2) = time_bounds(2)+3600;
    
    animation_time_series = linspace(0,max(time_bounds),frames);
    dt = animation_time_series(2) - animation_time_series(1);
    
    animation_time_series = sort(unique([animation_time_series, time_ending_series+1]));
    
    yaw_series = linspace(yaw_bounds(1),yaw_bounds(2)*yaw_revolutions,numel(animation_time_series));
    pitch_series = linspace(pitch_bounds(1),pitch_bounds(2)*pitch_revolutions,numel(animation_time_series));
    yaw_series = rem(yaw_series,360);
    pitch_series = rem(pitch_series,360);
    
    fade_tail_time = 3600*24;
    unit_size = 0.5e8;
    
    
    if record_video
        v = VideoWriter("MD_allocation_"+string(round(omega_end_range(omega_index)))+"_topview", 'MPEG-4');
        v.FrameRate = 60;
        open(v);
    end
    
    
    for ind_frame = 1 : numel(animation_time_series)
        fprintf("--------\n")
    
        t = animation_time_series(ind_frame);
        ind_frame
    
        scatter(nan,nan,"k")
        hold on
        grid on
        axis tight equal
    
        set(gca,"Position",[0,0.03,1,0.93])
    
        earth_r_check = [];
        luna_r_check = [];
    
        rocket_z_tracer = [];
        rocket_dot = [];
        rocket_tail = [];
        rocket_burn = [];
        lunar_ticks = [];
        rocket_label = [];
        extra = [];
    
        for n=1:numel(sim_labels)
            statematrix_spec = [];
            statematrix_spec = pso_results_total.(sim_labels(n)).results_struct.transit_statematrix;
    
            if t > min(statematrix_spec(:,1)) && t < max(statematrix_spec(:,1)) && ismember(pso_results_total.(sim_labels(n)).TA_omega(2), omega_end_range(omega_index))
    
                ind_transitions = find(diff(statematrix_spec(:,end)) ~= 0);
                statematrix_spec(ind_transitions+1,:) = [];
    
                [state_interp, ind_between] = interp_statematrix_via_time(statematrix_spec,t);
                [state_next_interp, ind_between_next] = interp_statematrix_via_time(statematrix_spec,t+dt);
    
                ind_range = [min([ind_between.', ind_between_next.']), max([ind_between.', ind_between_next.'])];
    
                burn_active = false;
                for m=1:numel(ind_transitions)
                    if ( ind_transitions(m) > ind_range(1) && ind_transitions(m) < ind_range(2) )
                        burn_active = true;
                        burn_v = state_next_interp(5:7) - state_interp(5:7);
                        break
                    end
                end
                if t+dt > statematrix_spec(end,1)
                    %final rendezvous burn
                    burn_active = true;
                    burn_v = state_next_interp(5:7);
                end
    
                if burn_active
                    rocket_dot(n) = scatter3(state_interp(2),state_interp(3),state_interp(4),10,"r","filled");
                else
                    rocket_dot(n) = scatter3(state_interp(2),state_interp(3),state_interp(4),10,"w","filled");
                end
                
                rocket_z_tracer(n) = plot3([state_interp(2),state_interp(2)],[state_interp(3),state_interp(3)],[0,state_interp(4)],"--w",LineWidth=0.01);
    
                if burn_active
                    burn_unit = burn_v./norm(burn_v);
                    rocket_burn(n) = quiver3(state_interp(2),state_interp(3),state_interp(4), burn_unit(1)*unit_size,burn_unit(2)*unit_size,burn_unit(3)*unit_size,"r","MaxHeadSize",0.9);
                else
                    rocket_burn = [];
                end
    
                %making fade tail
                tail_statematrix = [];
                [~,ind_tail_start] = min(abs(statematrix_spec(:,1)-(t-fade_tail_time)));
                [~,ind_tail_end] = min(abs(statematrix_spec(:,1)-t));
                tail_statematrix = statematrix_spec(ind_tail_start:ind_tail_end,:);
                tail_statematrix(end,:)=[];
                tail_statematrix = [tail_statematrix; state_interp];
                rocket_tail(n) = patch([tail_statematrix(:,2); nan],[tail_statematrix(:,3); nan],[tail_statematrix(:,4); nan],'w','EdgeColor','w',...
                    'FaceVertexAlphaData',linspace(0,0.5,height(tail_statematrix)+1).','AlphaDataMapping','none',...
                    'EdgeAlpha','interp');
    
                %where do these statematrixes say the earth and moon are?
                if interp1([min(statematrix_spec(:,1)), max(statematrix_spec(:,1))],[0,1],t) > 0.95
                    earth_r_check = [earth_r_check; state_interp(8:13)];
                    %luna_r_check = [luna_r_check; state_interp(14:19)];
                end
    
                TA_omega = pso_results_total.(sim_labels(n)).TA_omega;
                dv_total = pso_results_total.(sim_labels(n)).results_struct.dv_total;
                TW_ratio = pso_results_total.(sim_labels(n)).results_struct.acc_series;
                TW_ratio = max(TW_ratio)/9.81;
    
                rocket_desc = "MDR-"+ n + newline + "$\nu_l$ = " + TA_omega(1) + "$^{\circ}$" + newline;
                rocket_desc = rocket_desc + "$\Delta$v = "+ round(dv_total/1e3,3) + " km/s" + newline;
                rocket_desc = rocket_desc + "max T/W $\approx$ " + round(TW_ratio,3) + newline;
                %rocket_desc = "MDR-"+n;
    
                rocket_label(n) = text(state_interp(2),state_interp(3),state_interp(4), [" ", rocket_desc], Interpreter="latex", FontSize=10, VerticalAlignment="cap", Color = [1,1,1]);
            end
        end    
       
    
        if t < max(earth_init_statematrix(:,1))
            [earth_backup, ~] = interp_statematrix_via_time(earth_init_statematrix,t);
            earth_r_check = [earth_r_check; repelem( earth_backup(2:7),3,1) ]; %high weighting to init statematrix earth when its available
        end
        if t < max(lunar_init_statematrix(:,1)) 
            [luna_backup, ~] = interp_statematrix_via_time(lunar_init_statematrix,t);
            luna_r_check = [luna_r_check; repelem( luna_backup(2:7),3,1) ];
        end
    
        state_error=[];
        for n=1:width(earth_r_check)
            state_error(1,n) = max(earth_r_check(:,n)) - min(earth_r_check(:,n));
            state_error(2,n) = max(luna_r_check(:,n)) - min(luna_r_check(:,n));
            earth_state(n) = mean(earth_r_check(:,n));
            luna_state(n) = mean(luna_r_check(:,n));
        end
        luna_state(7) = rem(luna_backup(8),360);
    
        extra(1) = plot3(earth_init_statematrix(:,2),earth_init_statematrix(:,3),earth_init_statematrix(:,4),"b");
        extra(2) = plot3(lunar_init_statematrix(:,2),lunar_init_statematrix(:,3),lunar_init_statematrix(:,4),"w");
        
        [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(earth_state(1),earth_state(2),earth_state(3),earth_radius,earth_radius,earth_radius,64);
        extra(5) = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[0,0,1]);
        [surface_map_x,surface_map_y,surface_map_z] = ellipsoid(luna_state(1),luna_state(2),luna_state(3),lunar_radius,lunar_radius,lunar_radius,64);
        extra(6) = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[repelem(0.5,3)]);
    
        extra(3) = text(luna_state(1),luna_state(2),luna_state(3), [""," Luna" + newline + "$\nu = "+round(luna_state(7),2)+"^{\circ}$"], Interpreter="latex", FontSize=12,VerticalAlignment="cap",color=[1,1,1]);
        extra(4) = scatter3(luna_state(1),luna_state(2),luna_state(3),200,"w+");
    
        %plotting lunar orbit tick markers
        if view_type == 1
            lunar_tick_series = linspace(0,earthluna_orbital_period,36);
            lunar_tick_series(end) = [];
            for n=1:numel(lunar_tick_series)
                [lunar_orbit_tickstate, ~] = interp_statematrix_via_time(lunar_init_statematrix,lunar_tick_series(n));
                lunar_ticks(n) = plot3([lunar_orbit_tickstate(2),lunar_orbit_tickstate(2)],[lunar_orbit_tickstate(3),lunar_orbit_tickstate(3)],[0,lunar_orbit_tickstate(4)],"--",color=[repelem(0.5,3)],LineWidth=0.03);
            end
        end
    
    
        %plotting target orbit
        for n=1:numel(omega_end_range)
            if ismember(n,omega_index)
                target_inc = earth_obliquity + lunar_inc;
                target_omega = omega_end_range(n);
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
                ];
                target_r_initial = target_T_matrix*[0; target_radius; 0];
                target_v_initial = target_T_matrix*[-sqrt(G*earth_mass/target_radius); 0; 0];
        
                target_orbit_statematrix = [];
        
                target_orbit_stateinital = [target_r_initial.', target_v_initial.']; 
                target_orbit_timespan = [0,sqrt((target_radius^3)*(4*pi^2)/(G*earth_mass))];
                [target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),target_orbit_timespan,target_orbit_stateinital,odeset('Reltol',error_tolerance));
                target_orbit_statematrix = [target_timerange,target_state_matrix];
                for m=2:4
                    target_orbit_statematrix(:,m) = target_orbit_statematrix(:,m) + earth_state(m-1);
                end
                target_orbit(n) = plot3(target_orbit_statematrix(:,2),target_orbit_statematrix(:,3),target_orbit_statematrix(:,4),"w");
            end
        end

        xlim(x_bounds)
        ylim(y_bounds)
        zlim(z_bounds)
    
        ax = gca;
        ax.FontSize = 20;
    
        if view_type == 1
            view([yaw_series(ind_frame),pitch_series(ind_frame)])
            legend( "MD allocation phase transits $ \forall \enspace \Omega_t$ = "+ omega_end_range(omega_index)+"$^{\circ}$" + newline +"time: "+round(t/(3600*24),3) +" days", ...
                textcolor="w", FontSize=20, Interpreter="latex",location="northwest")
            legend boxoff
        else
            view([0,90]) 
        end
    
        set(gcf, 'Color', [1,1,1])
        set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')
    
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
    
        hold off
    
        drawnow()
    
        if record_video
            frame = getframe(gcf);
            for n=1:3
                frame_new(:,:,n) = uint8(imresize(squeeze(frame.cdata(:,:,n)), [nan, 1920*2 ]));
            end
            writeVideo(v,frame_new)
        end
   
        delete(rocket_dot);
        delete(rocket_tail);
        delete(target_orbit);
        delete(rocket_burn);
        delete(rocket_z_tracer);
        delete(target_orbit);
        delete(lunar_ticks);
        delete(rocket_label);
        delete(extra);
    end
    
    if record_video
        close(v);
        v = [];
    end

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