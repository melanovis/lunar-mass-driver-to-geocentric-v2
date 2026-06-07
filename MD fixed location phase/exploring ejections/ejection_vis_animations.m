format compact
clear
clc
clf reset

% ----------

record_video = true;

fps = 60;
frames = 200 * 60;

% time normalised
% lunar TA
% eject v
% azimuth
% elevation
% longitde
% latitude

keyframe_series = [
0, 0, 0, 90, 0, 0, 0
0.1, 0, 2.33e3, 90, 0, 0, 0
0.15, 0, 3e3, 90, 0, 0, 0
0.2, 0, 2.33e3, 90, 0, 0, 0
0.25, 1, 2.33e3, 90, 0, 0, 0
0.26, 1.25, 2.33e3, 90, 0, 0, 0
0.27, 1.25, 2.33e3, 0, 0, 0, 0
0.35, 1.25, 2.33e3, 360, 0, 0, 0
0.36, 1.25, 2.33e3, 90, 0, 0, 0
0.4, 1.25, 2.33e3, 90, 90, 0, 0
0.45, 2.25, 2.33e3, 0, 0, 0, 0
0.5, 2.25, 2.33e3, 90, 0, 0, 0
0.51, 2.25, 2.4e3, 90, 0, -45, 0
0.52, 2.25, 2.4e3, 90, 0, -45, 45
0.53, 2.25, 2.4e3, 90, 0, -45, -45
0.54, 2.25, 2.4e3, 90, 1e-3, -170, -60
0.55, 2.25, 2.4e3, 270, 1e-3, 170, 60
0.56, 2.25, 2.4e3, 90, 1e-3, -90, 0 %location revolving animation starts
0.65, 2.25, 2.4e3, 90, 1e-3, -90, 0 %location revolving animation starts
0.675, 2.25, 2.6e3, 90, 45, -90, 0
0.7, 2.25, 2.8e3, 270, 10, 90, 0
0.75, 1.25, 2.33e3, 45, 0.0148, -34.775, 9.378
0.8, 1.25, 3e3, 45, 0.0148, -34.775, 9.378
0.85, 2.25, 2.33e3, 45, 0.0148, -34.775, 9.378
0.9, 1.25, 2.33e3, 45, 0.0148, -34.775, 9.378
1, 0, 0, 90, 0, 0, 0
];
keyframe_series = single(keyframe_series);

[~,ind_sort] = sort(keyframe_series(:,1));
keyframe_series = keyframe_series(ind_sort,:);

keyframe_series(:,1) = keyframe_series(:,1)./max(keyframe_series(:,1));

pause_node_times = keyframe_series(:,1);

pause_frames = 5;

keyframe_duplicate = keyframe_series;
for n=1:pause_frames
    keyframe_series = [keyframe_series; [keyframe_duplicate(:,1)+n*(1/frames),keyframe_duplicate(:,2:end)]];
end

[~,ind_sort] = sort(keyframe_series(:,1));
keyframe_series = keyframe_series(ind_sort,:);

time_interp = linspace(0,1,frames);
for n=2:width(keyframe_series)
    keyframe_interp(:,n) = interp1(keyframe_series(:,1),keyframe_series(:,n),time_interp,"makima");
end
keyframe_interp(:,1) = time_interp;


%spinning animation
circ_range = [0.56, 0.65];
inds_circ = find( and(logical(time_interp > circ_range(1)), logical(time_interp < circ_range(2))) );
if mod(numel(inds_circ),2)==1
    inds_circ = [inds_circ, inds_circ(end)+1];
end
range_circ = [inds_circ(1),inds_circ(end)];
theta_pos_range = [linspace(0,360,numel(inds_circ)/2),linspace(0,360,numel(inds_circ)/2)];
theta_az_range = 90+[repelem(0,numel(inds_circ)/2),linspace(0,360*2,numel(inds_circ)/2)];
theta_az_range = rem(theta_az_range,360);
pos_start = [-90;0];
for n=1:numel(inds_circ)
    theta = theta_pos_range(n);
    pos_new(n,:) = [
    cosd(theta), sind(theta)
    -sind(theta), cosd(theta)
    ]*pos_start;
    az_new(n) = theta_az_range(n);
end
keyframe_interp(range_circ(1):range_circ(2),6:7) = pos_new;
keyframe_interp(range_circ(1):range_circ(2),4) = az_new;


pause_node_times = pause_node_times + (1/frames)*(pause_frames-1);
for n=1:numel(pause_node_times)
    [~,ind_pause] = min(abs(pause_node_times(n) - time_interp));
    pause_inds(n) = ind_pause;
end
pause_inds(1) = 1;




for n=2:numel(pause_inds)
    inds_range = [pause_inds(n-1):pause_inds(n)];
    fprintf("----------------------\n")

    if record_video
        v = [];
        v = VideoWriter("ejection_segment_"+string(n-1), 'MPEG-4');
        v.FrameRate = 60;
        open(v);
    end

    for m = 1:numel(inds_range)
        fprintf("-\n")
        interp_spec = double( keyframe_interp(inds_range(m),:) );
        interp_spec(1)
        interp_spec(2:end)
        %render_ejection(2, 2.5e3, 45, 0, -90, 0);
        render_ejection(interp_spec(2),interp_spec(3),interp_spec(4),interp_spec(5),interp_spec(6),interp_spec(7));
        if record_video
            frame = [];
            frame_new = uint8([]);
            frame = getframe(gcf);
            for n=1:3
                frame_new(:,:,n) = uint8(imresize(squeeze(frame.cdata(:,:,n)), [nan, 1920*2 ]));
            end
            writeVideo(v,frame_new)
        end
    end
    if record_video
        close(v);
    end
end


% ind_view = 4;
% hold on
% grid on
% plot(keyframe_interp(:,1),keyframe_interp(:,ind_view),"-k*")
% scatter(keyframe_series(:,1),keyframe_series(:,ind_view),"r","filled")
% scatter(keyframe_interp(pause_inds,1),keyframe_interp(pause_inds,ind_view),1000,"bx")