
%% #1 TrajGen.m 
% This document generates medial axes from a polygonal map
% It returns end point set "Ed", and end point coordinates EdPts
% and branch point set "Brch", and branch point coordinates BrchPts

%% #2 Planning.m
% This document generates a grid map on medial axes to obtain cost map from
% any locations to the goal, and returns a pixel-based cost map "Pathway" 

%% #3 EliminateFakeBranch.m
% This document eliminates end points that have been classfied as branch
% points wrongly, and returns two new sets "Ed0/EdPts0" and "Brch0/BrchPts0"

%% #4 TrajVec.m
% This generates paths from an end point/branch point to the next branch
% point, and records trajectories as connections of vectors. These vectors 
% have a maximum length of 1/2 channel width.

%% #5 trackplot.m
% plot robot trajectories with global control vector input 