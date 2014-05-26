% image acquisition
a=imread('101_1.tif');
% rotate image
a=imrotate(a,-20);
% write to a BitMap image the rotated mage
img=double(a);
% if image sizes are not 8*K where K is an integer
% image is resized

% call "centralizing" function
[out,xc,yc]=centralizing(img,0);
%[out,xc2,yc2]=centralizing2(img,0);

% show image and the center point determinated by
% the previous function
figure('Name','immagine');
imshow(a);
hold on;
plot(xc,yc,'O');
%plot(xc2,yc2,'X');
hold off;

