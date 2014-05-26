Unzip all files into Matlab current directory and type
"fprec" to start fingerprint image processing.
Type "helpwin fprec" for more detailed informations.


Notes:
old_version.zip is the previous release (Fingerprint Recognition System 2.0) 


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
