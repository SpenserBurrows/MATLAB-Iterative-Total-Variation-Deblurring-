function [u,i,cvg] = TVdeblur(z,k,a1,Threshold); %declare function
[mz1,nz1]=size(z); %get originial size of image
edgesize=round(max([mz1,nz1]/9)); %determine size of edge extrapolation

for i=1:edgesize;%extrapolate and taper off outer values to minimize edge effects
avg3=3+i*.0004;%taper off corners 
avg4=4+i*.0004;%taper off edges 
    
[m,n]=size(z); %get current size of z
ztl=(z(1,1)+z(2,1)+z(1,2))/avg3; %get top left corner
zbl=(z(end,1)+z(end-1,1)+z(end,2))/avg3; %get bottom left corner
zbr=(z(end,end)+z(end-1,end)+z(end,end-1))/avg3; %get bottom right corner
ztr=(z(1,end)+z(2,end)+z(1,end-1))/avg3;%get top right corner

zl=[ztl;(z(1:m-2,1)+z(2:m-1,1)+z(3:m,1)+z(2:m-1,2))/avg4;zbl];%get left edge
zr=[ztr;(z(1:m-2,n)+z(2:m-1,n)+z(3:m,n)+z(2:m-1,n-1))/avg4;zbr];%get right edge 
zt=[ztl,(z(1,1:n-2)+z(1,2:n-1)+z(1,3:n)+z(2,2:n-1))/avg4,ztr];%get top edge 
zb=[zbl,(z(m,1:n-2)+z(m,2:n-1)+z(m,3:n)+z(m-1,2:n-1))/avg4,zbr];%get bottom edge

z=[ztl,zt,ztr;zl,z,zr;zbl,zb,zbr];%form new image with edges 
end

mnstart=round((size(z)-[mz1,nz1])/2);%get beginning of original image inside new one
mnend=size(z)-mnstart-1;%get end of original image inside new one

g=1.5;%step size
ep=.001;%small constant to keep 1/gradient(u) term from blowing up
[mz,nz]=size(z); %get size of image
u=z;%use blurred image as initial guess
ki=flip(flip(k,1),2); %flip kernel
 
dfdu=zeros(mz,nz); %initialize arrays
TVterm=zeros(mz,nz);
convterm=zeros(mz,nz);
i=0;
ImageTitleFormat='Restored Image after %d iterations';
cvg=0;
avgdfdu=100;
 
%initialize image to show during iterations
    p1=subplot(1,2,1);
    image(z(mnstart(1)+1:mnend(1)+1,mnstart(2):mnend(2)));%show original
    title('Noisy Blurred Image');
 
    ImageTitle=sprintf(ImageTitleFormat,i);
    p2=subplot(1,2,2);
    image(u(mnstart(1)+1:mnend(1)+1,mnstart(2):mnend(2)));
    title(ImageTitle);
    truesize;
 
tic
while avgdfdu > Threshold %iterate until the mean value of df/du is below the threshold value
convterm=conv2((conv2(u,k,'same')-z),ki,'same'); %perform convolutions
TVterm=a1*del2(u)./sqrt(abs(gradient(u))+ep^2); %get TV term
dfdu=convterm-TVterm; %get df/du
u=u-g*dfdu; %get next u value
i=i+1; %increment i
avgdfdu=mean(mean(abs(dfdu)));%get mean value of df/du
cvg=[cvg,avgdfdu];%record mean value of df/du
 
ImageTitle=sprintf(ImageTitleFormat,i); %show image updating during iterations
p2=subplot(1,2,2);
image(u(mnstart(1)+1:mnend(1)+1,mnstart(2):mnend(2)));
title(ImageTitle);
drawnow;        
end
toc
cvg=cvg(2:end); %trim convergence series
u=u(mnstart(1)+1:mnend(1)+1,mnstart(2):mnend(2)); %trim output image back to original size

f2=figure;% show convergence 
plot([1:i],cvg,'k');
title('Convergance');
xlabel('number of iterations');
ylabel('Average value of df/du');
xlim([1,i])
end
