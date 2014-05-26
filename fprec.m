% Filterbank-Based Fingerprint Matching (A.K.Jain, S.Prabhakar, L.Hong and S.Pankanti, 2000)
%
% Abstract
% With identity fraud in our society reaching unprecedented proportions and with an increasing emphasis on 
% the emerging automatic personal identification applications, biometrics-based verification, especially 
% fingerprint-based identification, is receiving a lot of attention. There are two major shortcomings
% of the traditional approaches to fingerprint representation. For a considerable fraction of population, 
% the representations based on explicit detection of complete ridge structures in the fingerprint are 
% difficult to extract automatically. The widely used minutiae-based representation does not utilize a 
% significant component of the rich discriminatory information available in the fingerprints. Local ridge 
% structures cannot be completely characterized by minutiae. Further, minutiae-based matching has difficulty 
% in quickly matching two fingerprint images containing different number of unregistered minutiae points. 
% The proposed filter-based algorithm uses a bank of Gabor filters to capture both local and global details 
% in a fingerprint as a compact fixed length FingerCode. The fingerprint matching is based on the Euclidean 
% distance between the two corresponding FingerCodes and hence is extremely fast. We are able to achieve a 
% verification accuracy which is only marginally inferior to the best results of minutiae-based algorithms 
% published in the open literature. Our system performs better than a state-of-the-art minutiae-based system 
% when the performance requirement of the application system does not demand a very low false acceptance rate. 
% Finally, we show that the matching performance can be improved by combining the decisions of the matchers 
% based on complementary (minutiae-based and filter-based) fingerprint information. 
%
% Index Terms: Biometrics, FingerCode, fingerprints, flow pattern, Gabor filters, matching, texture, verification.
%
% Type "fprec" on Matlab command window to start image processing.
% This source code provides a new, improved GUI respect to the previous release.
% Simulation parameters can be changed in this file:
%
% n_bands: the number of concentric bands
% h_bands: the width of each band in pixels
% n_arcs: the number of arcs (each band has exactly n_arcs arcs)
% h_radius: the inner radius in pixel (the central band is not considered
%                                      because it is a too small area)
% num_disk: the number of gabor filters
%
% n_sectors and h_lato rapresents respectively the total number of sectors
% (the length of the feature vector associated to each filter of filter-bank:
% there are num_disk gabor filters) and the height of the cropped image in pixels.
% N_secors and h_lato should not be changed.
%
% 
%   $Revision: 1.0 $  $Date: 2002.10.02  $
%   $Revision: 2.0 $  $Date: 2003.11.29  $    
%   $Revision: 3.0 $  $Date: 2004.06.22  $   by Luigi Rosa
%                                            email:   luigi.rosa@tiscali.it 
%                                            mobile:  +393403463208
%                                            website: http://utenti.lycos.it/matlab
%
% 
%   Modified respect to the previous version:
%
%   - Major bugs fixed
%   - New GUI 
%   - 8 Gabor filters 0 22.5 45 67.5 90 112.5 135 157.5 degree
%   - Convolution is performed in frequency domain
%   - DataBase
%   - Fingerprint matching
%   - Error management
%   - Complex filtering techniques
%   - Improved core point determination
%   - Robustness against noise
%   - Modifiable simulation parameters 
%
% Input fingerprint should be 256 x 256 image 8-bit grayscale @ 500 dpi.
% If these conditions are not verified some parameters in m-functions
% should be changed in a proper way (such as, for example, Gabor filter
% parameters in gabor2d_sub function). See the cited references for more
% details.
%  
%   M-files included:
%   
%     -fprec.m:            this file. It initializes the entire image processing. The simulation
%                          parameters can be changed in this main file.
%     -centralizing.m:     a function which accept an input image and determines the coordinates
%                          of the core point. The core point is determinated by complex filtering.
%                          The region of interest is determinated fixing a minimum threshold value
%                          for the variance. Input image is divided into non-overlapping blocks and
%                          only blocks with a variance smaller than this threshold value are considered
%                          background. The logical matrix (associated to the region of interest) is first 
%                          closed (Matlab function imclose), then eroded (Matlab function imerode) with two
%                          given structuring elements. The image is "mirrored" before convolution with complex
%                          filter, then it is re-cropped to its original sizes.
%    -mirror.m:            a function which is used to "mirror" input image in  order to avoid undesired
%                          boundary effects (function used by centralizing.m).
%    -recrop.m:            a function used to resize the mirrored filtered image (function used by centralizing.m)
%    -conv2fft.m:          this function performs 2D FFT-based convolution. Type "help conv2fft" on Matlab command
%                          window for more details.
%    -whichsector.m:       a function used to determine (for each pixel of the cropped image) the corresponding
%                          sectors of the concentric bands (function used by sector_norm.m). 
%    -sector_norm.m:       a function used to normalize input image and to calculate the features vector
%    -cropping.m:          this function is used to cropp the input fingerprint image after the core point is
%                          determinated.
%    -gabor2d_sub.m:       a function used to calculate the coefficients of the gabor 2D filters.
%    -vedicentro.m:        this simple routines uses the M-function centralizing.m and it is used to display
%                          the core point.  
%
%
% A crucial step in fingerprint recognition is core point determination. 
% If any error occurs while cropping image you can use the auxiliary m-file
% "vedicentro.m": it visualizes the input fingerprint and the core point
% calculated by the m-function "centralizing.m".
%
%
%  Notes:
%  The computational load can be significantly reduced by recursive filtering techniques.
%  For a complete publication list of Lucas J. van Vliet please visit the following URL:
%  http://www.ph.tn.tudelft.nl/~lucas/publications/papersLJvV.html
%  Here you will find articles concernings a recursive implementation of the Gaussian filter,
%  of the derivative Gaussian filter and of Gabor filter.
% 
%  If you want to optimize the proposed method an excellent article is the following one:
%  Erian Bezhani, Dequn Sun, Jean-Luc Nagel, and Sergio Carrato, "Optimized filterbank 
%  fingerprint recognition", Proc. SPIE Intern. Symp. Electronic Imaging 2003, 20-24 
%  Jan. 2003, Santa Clara, California.  
%
%  This code was developed using:
%    MATLAB Version 6.5 and Image Processing Toolbox Version 3.2 (R13)
%    Operating System: Microsoft Windows 2000 Version 5.0 (Build 2195: Service Pack 4)
%    Java VM Version: Java 1.3.1_01 with Sun Microsystems Inc. Java HotSpot(TM) Client VM
%
%  Please contribute if you find this software useful.
%  Report bugs to luigi.rosa@tiscali.it
%   
%   
%   References:
%
%   Cheng Long Adam Wang, researcher
%   Fingerprint Recognition System
%   http://home.kimo.com.tw/carouse9/FRS.htm
%
%   A. K. Jain, S. Prabhakar, and S. Pankanti, "A Filterbank-based Representation for 
%   Classification and Matching of Fingerprints", International Joint Conference on 
%   Neural Networks (IJCNN), pp. 3284-3285, Washington DC, July 10-16, 1999. 
%   http://www.cse.msu.edu/~prabhaka/publications.html
%
%   "Fingerprint Classification and Matching Using a Filterbank", Salil Prabhakar
%   A DISSERTATION Submitted to Michigan State University in partial fulfillment 
%   of the requirements for the degree of DOCTOR OF PHILOSOPHY, Computer 
%   Science & Engineering, 2001
%   http://biometrics.cse.msu.edu/SalilThesis.pdf
%
%   Final Report 18-551 (Spring 1999) Fingerprint Recognition Group Number 19
%   Markus Adhiwiyogo, Samuel Chong, Joseph Huang, Weechoon Teo
%   http://www.ece.cmu.edu/~ee551/Old_projects/projects/s99_19/finalreport.html
%
%   Kenneth Nilsson and Josef Bigun, "Localization of corresponding points in 
%   fingerprints by complex filtering", Pattern Recognition Letters, 24 (2003) 2135-2144
%   School of Information Science, Computer and Electrical Engineering (IDE), Halmstad 
%   University, P.O. Box 823, SE-301 18, Halmstad, Sweden.
%
% 
% ************************************************************************
% This code required a lot of time to be developed. Please send the author
% money, food, drinks or improvements to the code itself.
% This is my postal address: 
% Luigi Rosa
% Via Centrale 35
% 67042 Civita di Bagno
% L'Aquila --- ITALY
%  
% mobile +39 340 3463208
% email luigi.rosa@tiscali.it
% website http://utenti.lycos.it/matlab
% *************************************************************************
%
%
%  


clear;
clc;
close all;
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

while chos~=possibility,
    chos=menu('Fingerprint Recognition System','Select image and add to database','Select image for fingerprint recognition','Info','Delete database',...
        'Fingerprint image: visualization','Gabor Filter: visualization','Exit');
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Calculate FingerCode and Add to Database
    if chos==1
        clc;
        close all;
        selezionato=0;
        while selezionato==0
            [namefile,pathname]=uigetfile({'*.bmp;*.tif;*.tiff;*.jpg;*.jpeg;*.gif','IMAGE Files (*.bmp,*.tif,*.tiff,*.jpg,*.jpeg,*.gif)'},'Chose GrayScale Image');
            if namefile~=0
                [img,map]=imread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('Select a grayscale image');
            end
            if (any(namefile~=0) && (~isgray(img)))
                disp('Select a grayscale image');
                selezionato=0;
            end
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
        
        message=strcat('FingerCode was succesfully added to database. Fingerprint no. ',num2str(fp_number));
        msgbox(message,'FingerCode DataBase','help');
    end
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Fingerprint recognition
    if chos==2
        clc;
        close all;
        selezionato=0;
        while selezionato==0
            [namefile,pathname]=uigetfile({'*.bmp;*.tif;*.tiff;*.jpg;*.jpeg;*.gif','IMAGE Files (*.bmp,*.tif,*.tiff,*.jpg,*.jpeg,*.gif)'},'Chose GrayScale Image');
            if namefile~=0
                [img,map]=imread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('Select a grayscale image');
            end
            if (any(namefile~=0) && (~isgray(img)))
                disp('Select a grayscale image');
                selezionato=0;
            end
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
        
        % memoria per feature vector d'ingresso
        vettore_in=zeros(num_disk*n_sectors,1);
        for (angle=0:1:num_disk-1)    
            gabor=gabor2d_sub(angle,num_disk);
            ComponentPrint=conv2fft(NormalizedPrint,gabor,'same');
            [disk,vector]=sector_norm(ComponentPrint,1);    
            finger_code{angle+1}=vector(1:n_sectors);
            vettore_in(angle*n_sectors+1:(angle+1)*n_sectors)=finger_code{angle+1};
        end
        
        
        
        % FingerCode of input fingerprint has just been calculated.
        % Checking with DataBase
        if (exist('fp_database.dat')==2)
            load('fp_database.dat','-mat');
            %---- alloco memoria -----------------------------------
            %...
            vettore_a=zeros(num_disk*n_sectors,1);
            vettore_b=zeros(num_disk*n_sectors,1);
            best_matching=zeros(fp_number,1);
            valori_rotazione=zeros(n_arcs,1);
            % start checking ---------------------------------------
            for scanning=1:fp_number
                fcode1=data{scanning,1};
                fcode2=data{scanning,2};
                for rotazione=0:(n_arcs-1)
                    p1=fcode1;
                    p2=fcode2;
                    % ruoto i valori dentro disco
                    for conta_disco=1:num_disk
                        disco1=p1{conta_disco};
                        disco2=p2{conta_disco};
                        for old_pos=1:n_arcs
                            new_pos=mod(old_pos+rotazione,n_arcs);
                            if new_pos==0
                                new_pos=n_arcs;
                            end
                            for conta_bande=0:1:(n_bands-1)
                                disco1r(new_pos+conta_bande*n_arcs)=disco1(old_pos+conta_bande*n_arcs);
                                disco2r(new_pos+conta_bande*n_arcs)=disco2(old_pos+conta_bande*n_arcs);
                            end
                        end
                        p1{conta_disco}=disco1r;
                        p2{conta_disco}=disco2r;
                    end
                    % ruoto i dischi circolarmente
                    for old_disk=1:num_disk
                        new_disk=mod(old_disk+rotazione,num_disk);
                        if new_disk==0
                            new_disk=num_disk;
                        end
                        pos=old_disk-1;
                        vettore_a(pos*n_sectors+1:(pos+1)*n_sectors)=p1{new_disk};
                        vettore_b(pos*n_sectors+1:(pos+1)*n_sectors)=p2{new_disk};                    
                    end
                    d1=norm(vettore_a-vettore_in);
                    d2=norm(vettore_b-vettore_in);
                    if d1<d2
                        val_minimo=d1;
                    else
                        val_minimo=d2;
                    end
                    valori_rotazione(rotazione+1)=val_minimo;
                end
                [minimo,posizione_minimo]=min(valori_rotazione);
                best_matching(scanning)=minimo;
            end
            [distanza_minima,posizione_minimo]=min(best_matching);
            beep;
            message=strcat('The nearest fingerprint present in DataBase which matchs input fingerprint is  : ',num2str(posizione_minimo),...
                ' with a distance of : ',num2str(distanza_minima));
            msgbox(message,'DataBase Info','help');
            
        else
            message='DataBase is empty. No check is possible.';
            msgbox(message,'FingerCode DataBase Error','warn');    
        end
        
    end % fine caso 2
    if chos==3
        clc;
        close all;
        helpwin fprec;
    end % fine caso 3
    if chos==4
        clc;
        close all;
        if (exist('fp_database.dat')==2)
            button = questdlg('Do you really want to remove the Database?');
            if strcmp(button,'Yes')
                delete('fp_database.dat');
                msgbox('Database was succesfully removed from the current directory.','Database removed','help');
            end
        else
            warndlg('Database is empty.',' Warning ')
        end
    end % fine caso 4
    if chos==5
        clc;
        close all;
        selezionato=0;
        while selezionato==0
            [namefile,pathname]=uigetfile({'*.bmp;*.tif;*.tiff;*.jpg;*.jpeg;*.gif','IMAGE Files (*.bmp,*.tif,*.tiff,*.jpg,*.jpeg,*.gif)'},'Chose GrayScale Image');
            if namefile~=0
                [img,map]=imread(strcat(pathname,namefile));
                selezionato=1;
            else
                disp('Select a grayscale image');
            end
            if (any(namefile~=0) && (~isgray(img)))
                disp('Select a grayscale image');
                selezionato=0;
            end
        end
        figure('Name','Selected image');
        imshow(img);
    end % fine caso 5
    if chos==6
        clc;
        close all;
        figure('Name','Gabor Filter');
        mesh(gabor2d_sub(0,num_disk));
    end % fine caso 6
end % fine while