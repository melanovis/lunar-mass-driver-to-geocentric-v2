format compact
clear
clc
clf reset

% ----------

G = 6.6743e-11;
error_tolerance = 1e-11;
earth_radius = 6.371e6; %m plus atmospheric margin
J_2 = 1.08262e-3;
earth_mass = 5.972e24;
mu = G*earth_mass;

orbit_inclination = 23.44 + 5.145;

parking_altitude = 3000e3;

target_radius = 1211203 + earth_radius;
parking_radius = parking_altitude + earth_radius;

target_orbital_period = sqrt(((target_radius^3)/mu)*(2*pi)^2);
parking_orbital_period = sqrt(((parking_radius^3)/mu)*(2*pi)^2);

target_precession = rad2deg( (-3/2) * ((earth_radius^2) / (target_radius^2)) * J_2 * (2*pi/target_orbital_period) * cosd(orbit_inclination) );
parking_precession = rad2deg( (-3/2) * ((earth_radius^2) / (parking_radius^2)) * J_2 * (2*pi/parking_orbital_period) * cosd(orbit_inclination) );

period_delta = abs(target_orbital_period - parking_orbital_period)
precession_delta = abs(target_precession - parking_precession);

max_wait_revs = target_orbital_period/period_delta;
max_wait_period = max_wait_revs*target_orbital_period

precession_omega_delta = precession_delta*max_wait_period;

parking_TA = 0;
target_TA = 180;

parking_T_matrix = [
1, 0, 0
0, cosd(orbit_inclination), -sind(orbit_inclination)
0, sind(orbit_inclination), cosd(orbit_inclination)
]*[
cosd(parking_TA), -sind(parking_TA), 0
sind(parking_TA), cosd(parking_TA), 0
0, 0, 1
];

depart_r = parking_T_matrix*[0; parking_radius; 0];
depart_v = parking_T_matrix*[-sqrt(G*earth_mass/parking_radius); 0; 0]; 

target_T_matrix = [
cosd(precession_omega_delta), -sind(precession_omega_delta), 0
sind(precession_omega_delta), cosd(precession_omega_delta), 0
0, 0, 1
]*[
1, 0, 0
0, cosd(orbit_inclination), -sind(orbit_inclination)
0, sind(orbit_inclination), cosd(orbit_inclination)
]*[
cosd(target_TA), -sind(target_TA), 0
sind(target_TA), cosd(target_TA), 0
0, 0, 1
];

arrive_r = target_T_matrix*[0; target_radius; 0];
arrive_v = target_T_matrix*[-sqrt(G*earth_mass/target_radius); 0; 0]; 


[v_1,v_2,delta_t] = lambert_solve_anticlockwise(depart_r, arrive_r, 1, mu);
transfer_dv = norm(v_1-depart_v) + norm(v_2-arrive_v)


transfer_orbit_timespan = [0, delta_t];
transfer_state_initial = [depart_r,v_1];
[transfer_timerange, transfer_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),transfer_orbit_timespan,transfer_state_initial,odeset('Reltol',error_tolerance));
transfer_orbit_statematrix = [transfer_timerange,transfer_state_matrix];

parking_orbit_timespan = [0, parking_orbital_period];
parking_state_initial = [depart_r,depart_v];
[parking_timerange, parking_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),parking_orbit_timespan,parking_state_initial,odeset('Reltol',error_tolerance));
parking_orbit_statematrix = [parking_timerange,parking_state_matrix];

target_orbit_timespan = [0, target_orbital_period];
target_state_initial = [arrive_r,arrive_v];
[target_timerange, target_state_matrix] = ode45(@(timerange, state_matrix)earth_twobody_ODE(timerange, state_matrix,earth_mass),target_orbit_timespan,target_state_initial,odeset('Reltol',error_tolerance));
target_orbit_statematrix = [target_timerange,target_state_matrix];



hold on
grid on
axis tight equal
[surface_map_x,surface_map_y,surface_map_z] = ellipsoid(0,0,0,earth_radius,earth_radius,earth_radius,64);
earth_obj = surf(surface_map_x,surface_map_y,surface_map_z,edgealpha=0.1,FaceColor=[0,0,1]);

plot3(parking_orbit_statematrix(:,2),parking_orbit_statematrix(:,3),parking_orbit_statematrix(:,4),"k")
plot3(target_orbit_statematrix(:,2),target_orbit_statematrix(:,3),target_orbit_statematrix(:,4),"k")
plot3(transfer_orbit_statematrix(:,2),transfer_orbit_statematrix(:,3),transfer_orbit_statematrix(:,4),"r")


view([30,30])


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