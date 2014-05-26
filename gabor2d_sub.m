function [gaborp_2d]=gabor2d_sub(angle,num_disk)
% Modified by Luigi Rosa
global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice


variance=32;
k=10;

x=cos(angle*pi/num_disk);
y=sin(angle*pi/num_disk);
w=(2*pi)/k;

p=0;
m=0;
for (i=-16:1:16)
   p=p+1;
   sinp(p)=i*y;
   cosp(p)=i*x;
   for (j=-16:1:16)
      m=m+1;
      x_s(m)=i;
      y_s(m)=j;
   end
end

p=0;
for (j=1:1:33)
   for (i=1:1:33)
      p=p+1;
      xx(p)=sinp(i)+cosp(j);
      yy(p)=cosp(i)-sinp(j);
      gaborp(p)=1*exp(-((xx(p)*xx(p))+(yy(p)*yy(p)))/variance)*cos(w*xx(p));
      gaborp_2d(i,j)=gaborp(p);
   end
end