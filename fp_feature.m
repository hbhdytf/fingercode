% clear;
% clc;
% close all;
% load('fp_database.dat','-mat');
% feature=data{1,1}{1,1}(:,1)';
% fd=fopen('feature.txt','wt');
% for m = 1:fp_number
% for n = 1:8
%     feature=data{m,1}{1,n}(:,1)';
%     for k = 1:80
%         if feature(k) >= 150
%             feature(k)=1;
%         else
%             feature(k)=0;
%         end
%         fprintf(fd,'%d',feature(k));
%     end
% end
% fprintf(fd,'\n');
% end
% fclose(fd);
clear;
clc;
close all;
%load('fp_database.dat','-mat');

%% For fprec
global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice num_disk


n_bands=5;
h_bands=20;
n_arcs=16;
h_radius=12;
h_lato=h_radius+(n_bands*h_bands*2)+16;
if mod(h_lato,2)==0
    h_lato=h_lato-1;
end
n_sectors=n_bands*n_arcs;
matrice=zeros(h_lato);
for ii=1:(h_lato*h_lato)
    matrice(ii)=whichsector(ii);
end
num_disk=8;
% 1--> add database
% 0--> recognition
%ok=0;
chos=0;
possibility=7;

messaggio='Insert the number of set: each set determins a class. This set should include a number of images for each person, with some variations in expression and in the lighting.';

pathname='H:\SHARE\Ground\DBII\DBII\';
fileFolder=fullfile(pathname);
dirOutput=dir(fullfile(fileFolder,'*.jpg'));
fileNames={dirOutput.name}';
fileN=length(fileNames);
delete('fp_database.dat');
for fm=1:fileN    
    namefile=char(fileNames(fm));
    if namefile~=0
        [img,map]=imread(strcat(pathname,namefile));
    end
    immagine=double(img);
    
    if isa(img,'uint8')
        graylevmax=2^8-1;
    end
    if isa(img,'uint16')
        graylevmax=2^16-1;
    end
    if isa(img,'uint32')
        graylevmax=2^32-1;
    end
    fingerprint = immagine;
    
    N=h_lato;
    
    [BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
    [CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
    [NormalizedPrint,vector]=sector_norm(CroppedPrint,0);
    
    for (angle=0:1:num_disk-1)
        gabor=gabor2d_sub(angle,num_disk);
        ComponentPrint=conv2fft(NormalizedPrint,gabor,'same');
        [disk,vector]=sector_norm(ComponentPrint,1);
        finger_code1{angle+1}=vector(1:n_sectors);
    end
    
    
    img=imrotate(img,180/(num_disk*2));
    fingerprint=double(img);
    
    [BinarizedPrint,XofCenter,YofCenter]=centralizing(fingerprint,0);
    [CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
    [NormalizedPrint,vector]=sector_norm(CroppedPrint,0);
    
    for (angle=0:1:num_disk-1)
        gabor=gabor2d_sub(angle,num_disk);
        ComponentPrint=conv2fft(NormalizedPrint,gabor,'same');
        [disk,vector]=sector_norm(ComponentPrint,1);
        finger_code2{angle+1}=vector(1:n_sectors);
    end
    % FingerCode added to database
    if (exist('fp_database.dat')==2)
        load('fp_database.dat','-mat');
        fp_number=fp_number+1;
        data{fp_number,1}=finger_code1;
        data{fp_number,2}=finger_code2;
        file{1,fp_number}=namefile;
        save('fp_database.dat','data','fp_number','file','-append');
    else
        fp_number=1;
        data{fp_number,1}=finger_code1;
        data{fp_number,2}=finger_code2;
        file{1,fp_number}=namefile;
        save('fp_database.dat','data','fp_number','file');
    end
    
    %message=strcat('FingerCode was succesfully added to database. Fingerprint no. ',num2str(fp_number));
    %msgbox(message,'FingerCode DataBase','help');
    %%
    S=regexp(namefile,'.jpg','split');
    str=char(S(1));
    outfeature=['f' str '.txt'];
    fpathname='H:\SHARE\Ground\DBII\FP\';
    fd=fopen(strcat(fpathname,outfeature),'wt');
    for n = 1:8
        feature=data{fp_number,1}{1,n}(:,1)';
        for k = 1:80
            if feature(k) >= 300
                a='11';
            elseif feature(k)>=200
                a='10';
            elseif feature(k)>=100
                a='01';
            else
                a='00';
            end
            fprintf(fd,'%s',a);
        end
    end
    fprintf(fd,'\n');
    
    fclose(fd);
end