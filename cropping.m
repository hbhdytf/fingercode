function [CroppedPrint] = cropping(XofCenter,YofCenter,CentralizedPrint)
% Modified by Luigi Rosa
global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice


N = h_lato;
M=size(CentralizedPrint,1);

imgN=size(CentralizedPrint,1);
imgM=size(CentralizedPrint,2);
%----------------------------------------------
% ELIMINATO
%  
%if (YofCenter+30) <= M
%   YofCenter = YofCenter + 20;
%else
%   YofCenter = M;
%end
%-----------------------------------------------
%center point can be shown on command line
%         ----------------
%         |              |
%         |              |
%         |       .      | y=row
%         |              |
%         |              |
%         ----------------
%                x=column
%XofCenter%   column of matrix
%YofCenter%   row of matrix
  
%-------------------------------------------------------------
%          if A= 1 2 3
%                4 5 6
%                7 8 9
%             B=A(1:2,2:3)
%              = 2 3
%                5 6
%       creates B by extracting the first twos and last two 
%       columns of A
%-------------------------------------------------------------
if (YofCenter-floor(N/2)<1)||(YofCenter+floor(N/2)>imgN)||(XofCenter-floor(N/2)<1)||(XofCenter+floor(N/2)>imgM)
    %message='Cropping error: when the input image is cropped an error occurs: a possible error during center point determination.';
    %msgbox(message,'Cropping Error','warn');
    % gestisco caso in cui il CroppedPrint esce fuori
    temp=zeros(imgN+2*h_lato,imgM+2*h_lato);
    temp(h_lato+1:h_lato+imgN,h_lato+1:h_lato+imgM)=CentralizedPrint;
    CroppedPrint=temp(YofCenter-floor(N/2)+h_lato:YofCenter+floor(N/2)+h_lato,XofCenter-floor(N/2)+h_lato:XofCenter+floor(N/2)+h_lato);
    return;
else
    CroppedPrint=CentralizedPrint(YofCenter-floor(N/2):YofCenter+floor(N/2),XofCenter-floor(N/2):XofCenter+floor(N/2));
    return;
end