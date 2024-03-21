clear all; close all; clc


%% Downsample factor - take every nth frame
nth = 2;

%% Convert videos
files(:,1) = uipickfiles('FilterSpec','*.mp4');
h = waitbar(0,['Processing... (1/', num2str(length(files)) ')']);

for jj = 1:length(files)
MP4name = ([files{jj}(1:end-4),'_',num2str(nth),'ds.mp4']);
v = VideoReader(files{jj});
total = v.NumFrames;
outputVideo = VideoWriter(MP4name,'MPEG-4');
outputVideo.FrameRate = round(v.FrameRate/nth,0);
open(outputVideo)

for kk = 1:nth:total
img = rgb2gray(read(v,kk));
writeVideo(outputVideo,img)
waitbar(kk/total,h,['Processing... (', num2str(jj), '/', num2str(length(files)) ')']);
end

close(outputVideo)

end

waitbar(1,h,'Conversions completed')
