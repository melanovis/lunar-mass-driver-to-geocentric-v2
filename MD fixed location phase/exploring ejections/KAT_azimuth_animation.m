format compact
clear
clc
clf reset

% ----------

record_video = true;

fps = 60;
frames = 70 * 60;

if record_video
    v = VideoWriter("KAT-3", 'MPEG-4');
    v.FrameRate = fps;
    open(v);
end

% time normalised
% lunar TA
% eject v
% azimuth
% elevation
% longitde
% latitude

% %kat 1
% keyframe_series = [
% 0, 0, 0, 90, 0, 0, 0
% 0.05, 0, 2.4e3, 90, 0.0148, -34.775, 9.378
% 0.25, 0, 2.4e3, 360+90, 0.0148, -34.775, 9.378
% 0.3, 0.3, 2.4e3, 360+90, 0.0148, -34.775, 9.378
% 0.5, 0.3, 2.4e3, 360*2+90, 0.0148, -34.775, 9.378
% 0.55, 0.5, 2.4e3, 360*2+90, 0.0148, -34.775, 9.378
% 0.75, 0.5, 2.4e3, 360*3+90, 0.0148, -34.775, 9.378
% 0.8, 0.7, 2.4e3, 360*3+90, 0.0148, -34.775, 9.378
% 0.95, 0.7, 2.4e3, 360*4+90, 0.0148, -34.775, 9.378
% 1, 1, 2.4e3, 360*4+90, 0, 0, 0
% ];

% %kat 2
% keyframe_series = [
% 0, 0, 0, 90, 0, 0, 0
% 0.05, 0.25, 2.33e3, 90, 0.0148, -34.775, 9.378
% 0.15, 0.25, 2.33e3, 90+360, 0.0148, -34.775, 9.378
% 0.2, 0.25, 2.4e3, 90+360, 0.0148, -34.775, 9.378
% 0.3, 0.25, 2.4e3, 90+360*2, 0.0148, -34.775, 9.378
% 0.35, 0.25, 2.5e3, 90+360*2, 0.0148, -34.775, 9.378
% 0.45, 0.25, 2.5e3, 90+360*3, 0.0148, -34.775, 9.378
% 0.5, 0.25, 2.7e3, 90+360*3, 0.0148, -34.775, 9.378
% 0.6, 0.25, 2.7e3, 90+360*4, 0.0148, -34.775, 9.378
% 0.65, 0.25, 2.9e3, 90+360*4, 0.0148, -34.775, 9.378
% 0.75, 0.25, 2.9e3, 90+360*5, 0.0148, -34.775, 9.378
% 0.8, 0.25, 3e3, 90+360*5, 0.0148, -34.775, 9.378
% 0.9, 0.25, 3e3, 90+360*6, 0.0148, -34.775, 9.378
% 1, 1, 0, 360*4+90, 0, 0, 0
% ];

%kat 3
keyframe_series = [
0, 0, 0, 90, 0, 0, 0
0.05, 0, 2.33e3, 45, 0.0148, -34.775, 9.378
0.15, 1, 2.33e3, 45, 0.0148, -34.775, 9.378
0.2, 1, 2.345e3, 45, 0.0148, -34.775, 9.378
0.3, 2, 2.345e3, 45, 0.0148, -34.775, 9.378
0.35, 2, 2.38e3, 45, 0.0148, -34.775, 9.378
0.45, 3, 2.38e3, 45, 0.0148, -34.775, 9.378
0.5, 3, 2.42e3, 45, 0.0148, -34.775, 9.378
0.6, 4, 2.42e3, 45, 0.0148, -34.775, 9.378
0.65, 4, 2.55e3, 45, 0.0148, -34.775, 9.378
0.75, 5, 2.55e3, 45, 0.0148, -34.775, 9.378
0.8, 5, 2.77e3, 45, 0.0148, -34.775, 9.378
0.9, 6, 2.77e3, 45, 0.0148, -34.775, 9.378
1, 7, 0, 90, 0, 0, 0
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


pause_node_times = pause_node_times + (1/frames)*(pause_frames-1);
for n=1:numel(pause_node_times)
    [~,ind_pause] = min(abs(pause_node_times(n) - time_interp));
    pause_inds(n) = ind_pause;
end
pause_inds(1) = 1;


for n=2:numel(pause_inds)
    inds_range = [pause_inds(n-1):pause_inds(n)];
    fprintf("----------------------\n")
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
end

if record_video
    close(v);
    sound(sin(2*pi*600*(0:1/14400:0.15)), 14400);
end


% ind_view = 4;
% hold on
% grid on
% plot(keyframe_interp(:,1),keyframe_interp(:,ind_view),"-k*")
% scatter(keyframe_series(:,1),keyframe_series(:,ind_view),"r","filled")
% scatter(keyframe_interp(pause_inds,1),keyframe_interp(pause_inds,ind_view),1000,"bx")