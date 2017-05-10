cd 'E:\SimuVideo2 T\SimuPics';


for ii = 1:n-1
   img_temp = uint8(cat(3,simuproc1(:,:,ii),simuproc2(:,:,ii),simuproc3(:,:,ii))); 
   filename = sprintf('TowRRTHRA%04d.png',ii);
   imwrite(img_temp,filename);
end