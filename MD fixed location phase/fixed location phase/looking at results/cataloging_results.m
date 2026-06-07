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

for n=1:height(sim_labels)
    struct_spec = pso_results_total.(sim_labels(n));
    TA_omega = struct_spec.TA_omega;
    dv_spec = struct_spec.results_struct.dv_total;
    acc_spec = max(struct_spec.results_struct.acc_series);
    eject_v = struct_spec.results_struct.MD_params(1);
    dt_total = struct_spec.results_struct.dt_total;
    results_list(n,:) = [round(TA_omega,2), dv_spec, acc_spec/9.81, eject_v, dt_total/(3600*24)]; 
end
results_list = single(results_list);

mean(results_list(:,3))
mean(results_list(:,4))
max(results_list(:,5))

omega_range = unique(results_list(:,2));
TA_range = unique(results_list(:,1));

for n=1:height(results_list)
    ind_ta = find(results_list(n,1) == TA_range);
    ind_o = find(results_list(n,2) == omega_range);
    dv_results_mat(ind_ta,ind_o) = results_list(n,3);
    TW_results_mat(ind_ta,ind_o) = results_list(n,4);
    eject_v_results_mat(ind_ta,ind_o) = results_list(n,5);
    time_results_mat(ind_ta,ind_o) = results_list(n,6);
end

cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
colormap(flip(cmap))

subplot(1,4,1)
hold on
grid on
axis tight equal
imagesc(eject_v_results_mat./1e3)
ax = gca;
ax.FontSize = 15;
h = colorbar;
set(get(h,'label'),'string','eject v (km/s)', Interpreter="latex", FontSize=20);
h.Color=[1,1,1];
setblack(ax)
yticks([1:numel(TA_range)])
xticks([1:numel(omega_range)])
xticklabels(cellstr(string(omega_range))+"^{\circ}")
yticklabels(string(round(TA_range)+"^{\circ}"))
xtickangle(90);
xlabel("target orbit longitude of ascending node $(\Omega)$", Interpreter="latex", FontSize=18)
ylabel("lunar true anomaly $(\nu)$", Interpreter="latex", FontSize=18)
threshold = mean([max(eject_v_results_mat(:)),min(eject_v_results_mat(:))]);
for n=1:height(eject_v_results_mat)
    for m=1:width(eject_v_results_mat)
        if eject_v_results_mat(n,m) < threshold
            colour_spec = [0,0,0];
        else
            colour_spec = [1,1,1];
        end
        text(m,n,string(round(eject_v_results_mat(n,m)/1e3,2)),color=colour_spec, HorizontalAlignment="center")
    end
end

subplot(1,4,2)
hold on
grid on
axis tight equal
imagesc(dv_results_mat./1e3)
ax = gca;
ax.FontSize = 15;
h = colorbar;
set(get(h,'label'),'string','rocket $\Delta$v (km/s)', Interpreter="latex", FontSize=20);
h.Color=[1,1,1];
setblack(ax)
yticks([1:numel(TA_range)])
xticks([1:numel(omega_range)])
xticklabels(cellstr(string(omega_range))+"^{\circ}")
yticklabels(string(round(TA_range)+"^{\circ}"))
xtickangle(90);
xlabel("target orbit longitude of ascending node $(\Omega)$", Interpreter="latex", FontSize=18)
ylabel("lunar true anomaly $(\nu)$", Interpreter="latex", FontSize=18)
threshold = mean([max(dv_results_mat(:)),min(dv_results_mat(:))]);
for n=1:height(eject_v_results_mat)
    for m=1:width(eject_v_results_mat)
        if dv_results_mat(n,m) < threshold
            colour_spec = [0,0,0];
        else
            colour_spec = [1,1,1];
        end
        text(m,n,string(round(dv_results_mat(n,m)/1e3,2)),Color=colour_spec, HorizontalAlignment="center")
    end
end


subplot(1,4,3)
hold on
grid on
axis tight equal
imagesc(TW_results_mat)
ax = gca;
ax.FontSize = 15;
setblack(ax)
h = colorbar;
set(get(h,'label'),'string','rocket max T/W (conservative)', Interpreter="latex", FontSize=20);
h.Color=[1,1,1];
setblack(ax)
yticks([1:numel(TA_range)])
xticks([1:numel(omega_range)])
xticklabels(cellstr(string(omega_range))+"^{\circ}")
yticklabels(string(round(TA_range)+"^{\circ}"))
xtickangle(90);
xlabel("target orbit longitude of ascending node $(\Omega)$", Interpreter="latex", FontSize=18)
ylabel("lunar true anomaly $(\nu)$", Interpreter="latex", FontSize=18)
threshold = mean([max(TW_results_mat(:)),min(TW_results_mat(:))]);
for n=1:height(eject_v_results_mat)
    for m=1:width(eject_v_results_mat)
        if TW_results_mat(n,m) < threshold
            colour_spec = [0,0,0];
        else
            colour_spec = [1,1,1];
        end
        text(m,n,string(round(TW_results_mat(n,m),2)),Color=colour_spec, HorizontalAlignment="center")
    end
end


subplot(1,4,4)
hold on
grid on
axis tight equal
imagesc(time_results_mat)
ax = gca;
ax.FontSize = 15;
setblack(ax)
h = colorbar;
set(get(h,'label'),'string','transit time (days)', Interpreter="latex", FontSize=20);
h.Color=[1,1,1];
setblack(ax)
yticks([1:numel(TA_range)])
xticks([1:numel(omega_range)])
xticklabels(cellstr(string(omega_range))+"^{\circ}")
yticklabels(string(round(TA_range)+"^{\circ}"))
xtickangle(90);
xlabel("target orbit longitude of ascending node $(\Omega)$", Interpreter="latex", FontSize=18)
ylabel("lunar true anomaly $(\nu)$", Interpreter="latex", FontSize=18)
threshold = mean([max(time_results_mat(:)),min(time_results_mat(:))]);
for n=1:height(eject_v_results_mat)
    for m=1:width(eject_v_results_mat)
        if time_results_mat(n,m) < threshold
            colour_spec = [0,0,0];
        else
            colour_spec = [1,1,1];
        end
        text(m,n,string(round(time_results_mat(n,m),2)),Color=colour_spec, HorizontalAlignment="center")
    end
end


md_desc = "fixed MD parameters: azimuth ($\theta$) = 45$^{\circ}$, elevation ($\epsilon$) = 0.0148$^{\circ}$, longitude ($\lambda$) = -34.775$^{\circ}$, latitude ($\phi$) = 9.378$^{\circ}$";

t = sgtitle(["luna-LEO MDR transit constellation results",md_desc],color=[1,1,1],fontsize=20, Interpreter="latex");

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')


% bins = 30;
% 
% subplot(1,3,1)
% hold on
% grid on
% histogram(results_list(:,5)./1e3,bins,FaceColor=[1,1,1],EdgeColor=[0,0,0],FaceAlpha=1);
% ax = gca;
% ax.FontSize = 20;
% setblack(ax)
% xlabel("MD eject v (km/s)", Interpreter="latex", FontSize=20)
% 
% subplot(1,3,2)
% hold on
% grid on
% histogram(results_list(:,3)./1e3,bins,FaceColor=[1,1,1],EdgeColor=[0,0,0],FaceAlpha=1);
% ax = gca;
% ax.FontSize = 20;
% setblack(ax)
% xlabel("rocket $\Delta$v (km/s)", Interpreter="latex", FontSize=20)
% 
% subplot(1,3,3)
% hold on
% grid on
% histogram(results_list(:,4),bins,FaceColor=[1,1,1],EdgeColor=[0,0,0],FaceAlpha=1);
% ax = gca;
% ax.FontSize = 20;
% setblack(ax)
% xlabel("rocket max T/W (conservative)", Interpreter="latex", FontSize=20)
% 
% set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')


im_raw = getframe(gcf);
imwrite(im_raw.cdata, "figure.png"); 

function setblack(ax)
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
end