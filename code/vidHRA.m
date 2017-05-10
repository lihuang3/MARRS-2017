workdir = 'E:\SimuVideo2 T\SimuPics';
cd(workdir)
fig_list = dir('*png');
fps = 60;
filename = sprintf('vid_fps%03d.avi',fps);
outputVideo = VideoWriter(filename);
outputVideo.FrameRate = fps;

open(outputVideo)
temp_img =imread(fig_list(1).name);
ImCol = size(temp_img,2);
ImRow = size(temp_img,1);
numEle = numel(temp_img);

for ii = 1:length(fig_list)
%    if mod(ii,100)==0
      img = imread(fig_list(ii).name);
      writeVideo(outputVideo,img)
%    end
end

close(outputVideo)
