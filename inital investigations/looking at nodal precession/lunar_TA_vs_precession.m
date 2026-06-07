format compact
clear
clc
clf reset

% ----------

run_sweep = false;

G = 6.6743e-11;
earth_radius = 6.371e6;
earth_mass = 5.972e24;

lunar_precession_period = 6798.383*3600*24; %secs

mu_single = G*earth_mass;

station_orbit_inc = single( 23.44 + 5.145 );

res = 3e3;
alt_range = linspace(1200e3,1300e3,res);
ecc_range = linspace(0,0.5,res);

if run_sweep
    for ind_a = 1:numel(alt_range)
        for ind_e = 1:numel(ecc_range)
            
            station_orbit_radius = single( alt_range(ind_a) + earth_radius ); %periapse
            station_orbit_ecc = ecc_range(ind_e);
            station_orbit_semimajor = station_orbit_radius / (1-station_orbit_ecc);
            station_orbit_period = 2*pi*sqrt((station_orbit_semimajor^3)/mu_single);
            
            J_2 = 1.08262e-3;
            station_precession_rate = (-3/2) * ((earth_radius^2) / ((station_orbit_semimajor*(1-station_orbit_ecc^2))^2)) * J_2 * (2*pi/station_orbit_period) * cosd(station_orbit_inc);
            station_precession_rate = rad2deg(station_precession_rate);
            
            precession_cycles = floor(abs(station_precession_rate)*lunar_precession_period / 360);
            precession_remainder = rem(abs(station_precession_rate)*lunar_precession_period,360);
            
            rem_series(ind_e,ind_a) = precession_remainder;
            cycles_series(ind_e, ind_a) = precession_cycles;

        end
        if rem(ind_a,100)==0
            fprintf("- %i.\n",ind_a)
        end
    end

    revs_whole = [];

    for m=1:height(rem_series)
        for n=2:width(rem_series)-1
            if rem_series(m,n) < rem_series(m,n+1)
                revs_whole = [revs_whole ; alt_range(n),ecc_range(m),rem_series(m,n),cycles_series(m,n)];
            end
        end
    end

    save("nodal_results.mat","rem_series","cycles_series","alt_range","ecc_range","revs_whole")

end
load("nodal_results.mat")


cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
colormap(cmap)

hold on
grid on
axis tight square
scatter(revs_whole(:,1)./1e3,revs_whole(:,2),5,"w","filled")
ax = gca;
ax.FontSize = 20;
xlim([min(alt_range),max(alt_range)]./1e3)

xlabel("target altitude periapse (km)", Interpreter="latex", FontSize=20)
ylabel("target eccentricty", Interpreter="latex", FontSize=20)

title("SLSO orbits in the 1200-1300km altitude range", Interpreter="latex", FontSize=20)

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

%xlabel("target altitude (km)", Interpreter="latex", FontSize=20)
