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

parking_altitude_range = linspace(500e3,2000e3,2e3);

target_altitude = 1211203;
target_radius = target_altitude + earth_radius;

for ind_p = 1:numel(parking_altitude_range)

    parking_altitude = parking_altitude_range(ind_p);
    
    parking_radius = parking_altitude + earth_radius;
    
    target_orbital_period = sqrt(((target_radius^3)/mu)*(2*pi)^2);
    parking_orbital_period = sqrt(((parking_radius^3)/mu)*(2*pi)^2);
    
    target_precession = rad2deg( (-3/2) * ((earth_radius^2) / (target_radius^2)) * J_2 * (2*pi/target_orbital_period) * cosd(orbit_inclination) );
    parking_precession = rad2deg( (-3/2) * ((earth_radius^2) / (parking_radius^2)) * J_2 * (2*pi/parking_orbital_period) * cosd(orbit_inclination) );
    
    period_delta = abs(target_orbital_period - parking_orbital_period);
    precession_delta = abs(target_precession - parking_precession);
    
    max_wait_revs = target_orbital_period/period_delta;
    max_wait_period = max_wait_revs*target_orbital_period;
    
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
    
    [v_1,v_2,delta_t] = lambert_solve_anticlockwise(depart_r, arrive_r, 1, mu); %ideal dv for circular orbits
    transfer_dv = norm(v_1-depart_v) + norm(v_2-arrive_v);

    max_wait_profile(ind_p) = max_wait_period;
    dv_profile(ind_p) = transfer_dv;
end

[~,ind_mark] = mink(abs(dv_profile - 100),2);
parking_altitude_range(ind_mark)
max_wait_profile(ind_mark)./3600
dv_profile(ind_mark)

subplot(1,2,1)
hold on
grid on
axis tight
plot([target_altitude,target_altitude]./1e3,[min(max_wait_profile),max(max_wait_profile)]./3600,"--w",LineWidth=1.5)

plot(parking_altitude_range./1e3,max_wait_profile./3600,"w",LineWidth=1.5)
ax = gca;
ax.FontSize = 20;
set(gca,"YScale","log")
xlabel("parking altitude (km)", Interpreter="latex", FontSize=20)
ylabel("max rendezvous waiting period (hours)", Interpreter="latex", FontSize=20)
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

legend("target orbit altitude", Interpreter="latex", FontSize=18,TextColor=[1,1,1])
legend boxoff


subplot(1,2,2)
hold on
grid on
axis tight
plot(parking_altitude_range./1e3,dv_profile,"w",LineWidth=1.5)
plot([target_altitude,target_altitude]./1e3,[min(dv_profile),max(dv_profile)],"--w",LineWidth=1.5)
ax = gca;
ax.FontSize = 20;
xlabel("parking altitude (km)", Interpreter="latex", FontSize=20)
ylabel("max transfer-rendezvous $\Delta$v (m/s)", Interpreter="latex", FontSize=20)
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

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

im_raw = getframe(gcf);
imwrite(im_raw.cdata, "test.jpg"); 
