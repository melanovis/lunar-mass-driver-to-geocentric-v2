format compact
clear
clc
clf reset

%----

lunar_obliquity = 6.68;

load("pso_results_collected.mat")
load("craters_formatted.mat")

sim_labels = string(fieldnames(pso_results_total));

for n=1:numel(sim_labels)
    dv_series(n) =  pso_results_total.(sim_labels(n)).results_struct.dv_total;
    acc_series(n) = max(pso_results_total.(sim_labels(n)).results_struct.acc_series);
end
dv_series = sort(dv_series);

hold on
grid on
axis tight equal
xlim([-180,180])
ylim([-90,90])

surfacemap = imread('lunar_map.png');
surfacemap = flipud(surfacemap);

wm = image(surfacemap,'xdata',[-180,180],'ydata',[-90 90]);

set(gca,'ydir','normal')
uistack(wm,'down')

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
cmap = flip(cmap);

unit_vec_size = 6;

for n=numel(sim_labels):-1:1
    
    struct_spec = pso_results_total.(sim_labels(n));
    MD_params_spec = struct_spec.results_struct.MD_params;

    target_params = struct_spec.TA_omega;
    
    spec_lat = MD_params_spec(3);
    spec_long = MD_params_spec(2);
    spec_azimuth = MD_params_spec(4);
    spec_elevation = MD_params_spec(5);
    dv_spec = struct_spec.results_struct.dv_total;
    acc_spec = max(struct_spec.results_struct.acc_series);

    ind_colour = round(interp1([min(dv_series),max(dv_series)],[1,height(cmap)],dv_spec));
    %ind_colour = round(interp1([min(acc_series),max(acc_series)],[1,height(cmap)],acc_spec));

    eject_vector = [
    cosd(spec_azimuth), sind(spec_azimuth)
    -sind(spec_azimuth), cosd(spec_azimuth)
    ]*[0;1];

    azimuth_series(n) = spec_azimuth;
    long_series(n) = spec_long;
    lat_series(n) = spec_lat;
    elevation_series(n) = spec_elevation;

    quiver(spec_long,spec_lat,eject_vector(1)*unit_vec_size,eject_vector(2)*unit_vec_size,"w",MaxHeadSize=1)
    scatter(spec_long,spec_lat,150,"wx")
    scatter(spec_long,spec_lat,"filled",markerfacecolor = cmap(ind_colour,:), markeredgecolor = [0,0,0])

end

centroid = [mean(long_series), mean(lat_series)];

for n=1:numel(long_series)
    mean_sep(n) = norm(centroid - [long_series(n),lat_series(n)]);
end
std_sep = std(mean_sep);

% std_region = rectangle('Position', [centroid-std_sep/2,std_sep,std_sep], 'Curvature', [1 1], facecolor = [1,1,1,0.3],Edgecolor = [1,1,1,0.5]);
% scatter(centroid(1),centroid(2),"r+")


% for n=1:height(crater_labels)
%     if crater_latlong(n,2) < 170 && norm([crater_latlong(n,:) - flip(centroid)]) < 50 && crater_diameter(n) > 20
%         scatter(crater_latlong(n,2),crater_latlong(n,1),"rx")
%         text(crater_latlong(n,2),crater_latlong(n,1)," "+string(crater_labels(n)),fontsize = 9,Color=[1,1,1])
%     end
% end

ax = gca;
ax.FontSize = 20;

% xlim([-50,50]+centroid(1))
% ylim([-25,25]+centroid(2))

colormap(cmap)

h = colorbar;
set(get(h,'label'),'string','transit rocket $\Delta$v (km/s)', Interpreter="latex", FontSize=20);
%set(get(h,'label'),'string','conservative transit rocket acceleration m/s$^2$', Interpreter="latex", FontSize=20);
h.Color = [1,1,1];
clim([min(dv_series),max(dv_series)]./1e3)
%clim([min(acc_series),max(acc_series)])

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
imwrite(im_raw.cdata, "figure.jpg"); 