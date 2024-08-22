function LVMLTVDB(filename,color,blurType,Threshold,length,angle); %declare function
[Original]=imread(filename); %read in image
[m,n,c]=size(Original); %get image size
if (c==3) %if a color image
[Original,OrigMap]=imread(filename); %get colormap
else
GrayOriginal=Original; %if grayscale, use orignal image
end

a1=.001 %set a1 value

if (blurType == 1) %if out of focus blur
k=fspecial('disk',length); %define out of focus blurring kernel
else
k=fspecial('motion',length,angle); %define motion blurring kernel
end

if (color == 0) | (c==1) %if grayscale selected or if the image is greyscale to begin with
    if c==3
GrayOriginal=(rgb2gray(Original)); %convert to grayscale if original is RGB
    end
figure('Name','Grayscale');%open figure
colormap(gray(256));%set colormap
u=TVdeblur(double(GrayOriginal),k,a1,Threshold);%execute deblurring
newfile=strcat(erase(filename,'.jpg'),'_restored.jpg')%generate restored image filename
imwrite(u,gray(256),newfile);%save restored image

else %if RGB

[urgb,cvr,cvg,cvb]=TVdeblurRGB(Original,k,a1,Threshold);%execute RGB deblurring

ired=cvr(1);%get iterations for red, green, and blue
igreen=cvg(1);
iblue=cvb(1);
cvgred=cvr(2:end);%get convergence sequences for red, green, and blue
cvggreen=cvg(2:end);
cvgblue=cvb(2:end);
f2=gcf;%get tag for current figure
f2.Name='RGB';%plot final blurred and deblurred images
    p1=subplot(1,2,1);
    image(uint8(Original));
    title('Noisy Blurred Image');
colormap(OrigMap);
    itotal=ired+igreen+iblue;
    ImageTitleFormat='Restored Image after %d iterations';
    ImageTitle=sprintf(ImageTitleFormat,itotal);
    p2=subplot(1,2,2);
    image(uint8(urgb));
    title(ImageTitle);
colormap(OrigMap);
truesize;


figure('Name','Convergence','Position',[200,50,800,200]);%create convergence figure

p3=subplot(1,3,1);%plot red convergence
plot([1:ired],cvgred,'r');
title('Red Convergence');
xlabel('number of iterations');
ylabel('Average value of df/du');
xlim([1,ired]);

p4=subplot(1,3,2);%plot green convergence
plot([1:igreen],cvggreen,'g');
title('Green Convergence');
xlabel('number of iterations');
ylabel('Average value of df/du');
xlim([1,igreen]);

p5=subplot(1,3,3);%plot blue convergence
plot([1:iblue],cvgblue,'b');
title('Blue Convergance');
xlabel('number of iterations');
ylabel('Average value of df/du');
xlim([1,iblue]);

newfile=strcat(erase(filename,'.jpg'),'_restored.jpg')%generate restored image filename
imwrite(uint8(urgb),colormap(OrigMap),newfile);%save restored image

end
end
