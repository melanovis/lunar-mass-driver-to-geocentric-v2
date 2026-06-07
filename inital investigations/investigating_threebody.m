format compact
clear
clc
clf reset

% ----------

G = 6.6743e-11;
error_tolerance = 1e-10;

lunar_radius = 1.7374e6;
earth_radius = 6.371e6; %m

lunar_mass = 7.35e22; %kg
earth_mass = 5.972e24;
mu = G*(earth_mass + lunar_mass);

approx_semimajor = 384396e3;

pair_eccentricity = 0.0549;

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

timespan = linspace(0,earthluna_orbital_period*2,3e3);

[timerange, statematrix] = ode45(@(timerange, statematrix) earth_moon_ODEs(timerange, statematrix, earth_mass, lunar_mass), timespan, earthlunar_state_intial ,odeset('Reltol',error_tolerance));
earth_statematrix = [timerange, statematrix(:,1:6)];
lunar_statematrix = [timerange, statematrix(:,7:12)];

lunar_peri_state = lunar_statematrix(1,:);
earth_apo_state = earth_statematrix(1,:);

lunar_TA_list = get_planet_TAlist(lunar_peri_state(2:4),lunar_peri_state(5:7),lunar_statematrix);
earth_TA_list = get_planet_TAlist(earth_apo_state(2:4),earth_apo_state(5:7),earth_statematrix);
earth_TA_list = rem(earth_TA_list+180,360); %TA function is designed to start at periapse, accounting for earth's apo start

earth_TA_list = accumulate_increasing_TAlist(earth_TA_list);
lunar_TA_list = accumulate_increasing_TAlist(lunar_TA_list);

earth_statematrix = [earth_statematrix, earth_TA_list.'];
lunar_statematrix = [lunar_statematrix, lunar_TA_list.'];



subplot(1,2,1)
hold on
grid on
axis equal
scatter3(earth_r_initial(1),earth_r_initial(2),earth_r_initial(3),5,"b","filled")
scatter3(lunar_r_initial(1),lunar_r_initial(2),lunar_r_initial(3),5,"k","filled")

unit_vec_size = 5e7;
lunar_v_norm = lunar_v_initial./norm(lunar_v_initial);
earth_v_norm = earth_v_initial./norm(earth_v_initial);

quiver3(earth_r_initial(1),earth_r_initial(2),earth_r_initial(3), earth_v_norm(1)*unit_vec_size, earth_v_norm(2)*unit_vec_size, earth_v_norm(3)*unit_vec_size, "r")
quiver3(lunar_r_initial(1),lunar_r_initial(2),lunar_r_initial(3), lunar_v_norm(1)*unit_vec_size, lunar_v_norm(2)*unit_vec_size, lunar_v_norm(3)*unit_vec_size, "r")

plot3(earth_statematrix(:,2),earth_statematrix(:,3),earth_statematrix(:,4),"b")
plot3(lunar_statematrix(:,2),lunar_statematrix(:,3),lunar_statematrix(:,4),"k")

view([30,30])


subplot(1,2,2)
for n=1:50:height(lunar_statematrix)

    scatter(nan,nan)
    hold on
    grid on
    axis equal
    scatter3(earth_statematrix(n,2),earth_statematrix(n,3),earth_statematrix(n,4),10,"b","filled")
    scatter3(lunar_statematrix(n,2),lunar_statematrix(n,3),lunar_statematrix(n,4),10,"k","filled")
    plot3(earth_statematrix(:,2),earth_statematrix(:,3),earth_statematrix(:,4),"b")
    plot3(lunar_statematrix(:,2),lunar_statematrix(:,3),lunar_statematrix(:,4),"k")

    text(earth_statematrix(n,2),earth_statematrix(n,3),earth_statematrix(n,4)," "+string(earth_statematrix(n,8)))
    text(lunar_statematrix(n,2),lunar_statematrix(n,3),lunar_statematrix(n,4)," "+string(lunar_statematrix(n,8)))

    view([0,90])
    drawnow()

    hold off
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

function state_out = earth_moon_ODEs(t, state_in, earth_mass, lunar_mass)

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