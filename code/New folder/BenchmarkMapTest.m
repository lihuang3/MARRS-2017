clear;
close all;
clc

mapname = {'figT','figStdMap','figVB', 'maze','maze1'};
scale = {1, 0.5, 1, 1, 1};
%goalloc = {[245,273],[331,333],[297,169]};

goalloc = {[280,327],[331,333],[297,169],[615, 620],[490 308]};
dist_threshold = {[35 35],[35 35], [15 15],[15 15], [12 12]};
max_step = {15, 8, 10, 10, 8};
mission_increment = {1/2, 1/2, 1/3, 1/3, 1/3};
maps = cell(5,1);
for ii = 1:5
    maps{ii} = struct('name', mapname(ii), 'goal_loc', goalloc(ii), ...
        'distance_threshold', dist_threshold(ii),'magnification', ...
        scale(ii),'max_step',max_step(ii),'mission_increment',mission_increment(ii));
end

addpath(genpath(pwd));

fprintf('Choose a map number from the list ...\n')
for ii = 1:size(mapname,2)
   fprintf('Map #%d: %s\n', ii, mapname{ii});
end

temp = input('Input the map name number: ');
testmap = strcat(maps{temp,1}.name,'.fig');



while ~exist(testmap)
    fprintf('File name does not exist!\n')
    fprintf('Choose a map number from the list ...\n')
    for ii = 1:size(mapname,2)
       fprintf('Map #%d: %s\n', ii, mapname{ii});
    end
    temp = input('Input the map name number: ');
    testmap = strcat(maps{temp,1}.name,'.fig');
end

fprintf('Map %s is chosen.\n', testmap);
output = TestObj(maps{temp});

% filename = strcat(maps{temp,1}.name,'_output.mat');
% fprintf('Writing data to %s ... \n',filename)
% save(filename,'output')
% fprintf('Data saved!\n')



