close all
clear
clc

%cd ('E:\SimuVideo2 T\Alg1')
cd('C:\Users\lhuang28\Desktop\SimuVideo2 T\Alg1')
%load('dataT_owRRT1212.mat');
open('figStdMap.fig');

fig = getframe;
bw = imbinarize(fig.cdata(:,:,1),0.9);
BW = imresize(bw,0.5);
[l10,l20] = size(BW);
l1 = uint16(l10*1.2);
l2 = uint16(l20*1.2);
BW0 = uint16(zeros(l1,l2));
BW0(uint16((l1-l10)/2):uint16((l1-l10)/2)-1+l10,uint16((l2-l20)/2):l20-1+uint16((l2-l20)/2))=...
    BW;
BW = logical(BW0);
EBW = ~edge(BW,'canny');
Skel = bwmorph(BW,'skel',Inf);
Brch = bwmorph(Skel,'branchpoints');
[row col] = find(Brch);
BrchPts = [row col]; 
Ed = bwmorph(Skel,'endpoints');
[row col] = find(Ed);
EdPts = [row col];
figure
imshow(Skel);
hold on
%plot(BrchPts(:,2),BrchPts(:,1),'r.');
f1 = scatter(BrchPts(:,2),BrchPts(:,1),20,'filled');
f1.CData = [1 0 0];
f2 = scatter(EdPts(:,2),EdPts(:,1),20,'filled');
f2.CData = [0 0 1];
%plot(EdPts(:,2),EdPts(:,1),'b.');



