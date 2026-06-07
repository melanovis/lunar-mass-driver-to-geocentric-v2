format compact
clear
clc
clf reset

%-------

G = 6.6743e-11;

earth_radius = 6.371e6;
earth_mass = 5.972e24;

orbit_radius = 21054.35e3 + earth_radius; %circular orbit
orbit_inclination = 23.44 + 5.145;

mu = G*earth_mass;
orbit_period = 2*pi*sqrt((orbit_radius^3)/mu);

lunar_precession_period = 6798.383*3600*24; %secs

J_2 = 1.08262e-3;

precession_rate = (-3/2) * ((earth_radius^2) / (orbit_radius^2)) * J_2 * (2*pi/orbit_period) * cosd(orbit_inclination);

precession_rate_deg = rad2deg(precession_rate);
precession_daily = precession_rate_deg*3600*24;

abs(precession_rate_deg)*lunar_precession_period
