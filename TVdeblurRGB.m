function [u,cvr,cvg,cvb] = TVdeblurRGB(z,k,a1,Threshold); %declare function
[mz,nz,c]=size(z);%get size of image
z=double(z);%convert image to double representation 
[zR,zG,zB]=imsplit(z);%split color channels
u=zeros(mz,nz,c);%initalize output image

colormap(1/255*cat(2,[1:255]',zeros(2,255)'));%do red
[u(1:end,1:end,1),ired,cvgred]=TVdeblur(zR,k,a1,Threshold);
cvr=[ired,cvgred];

colormap(1/255*cat(2,zeros(1,255)',[1:255]',zeros(1,255)'));%do green
[u(1:end,1:end,2),igreen,cvggreen]=TVdeblur(zG,k,a1,Threshold);
cvg=[igreen,cvggreen];

colormap(1/255*cat(2,zeros(2,255)',[1:255]'));%do blue
[u(1:end,1:end,3),iblue,cvgblue]=TVdeblur(zB,k,a1,Threshold);
cvb=[iblue,cvgblue];
end