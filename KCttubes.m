%%  Basic Sub-resolution TT Diameter Method
%   Description: Example processing pipeline for estimating sub-resolution 
%   t-tubule diameters from confocal images as described in 
%   Kong, C.H.T. & Cannell, M.B. (2026) J. Microsc., and first demonstrated
%   in Soeller, C. & Cannell, M.B. (1999) Circ. Res. 84, 266–275.
%   Code tested in Matlab 2023b
%
%   Note that this is an example processing strategy. Exact algorithms 
%   and steps should be checked and modified to suit user data. 
%
%   Copyright (c) 2026, M.B. Cannell, C.H.T. Kong.
%   Permission is hereby granted, free of charge, to any person obtaining
%   a copy of this software and associated documentation files (the
%   "Software"), to use the Software without restriction, including
%   the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, but not sell copies of the Software.

%%  1.  Clear workspace and set options
clear all, close all                            % (optional) clear workspace
fwhmXY=0.35;                                    % PSF FWHM in xy, micron
fwhmZ=0.7;                                      % PSF FWHM in z, micron
colmap=colormap("jet");                         % colormap for visualisation

%%  2.  Read in some low S/N tt data and print key properties
load TTAnalysis_ExampleData05.mat
datSz=size(data);
fprintf("Data size: y=%i; x=%i; z=%i px.\n",datSz);
fprintf("Voxel size: xy=%1.2f; z=%1.2f microns.\n",dxy,dz);

%%  3.  Pre-processing 
%   3a. Data: Remove instrument dc offset (set by user). Do any other 
%   needed data pre-processing here such as cropping etc.
data=data-100;

%%  4.  PSF and masks for cell
%   4a. Generate PSF for deconvolution: Here, the PSF is approximated by a 
%   3D Gaussian function, however, more sophisticated calculations can be
%   used, e.g. the PSF Generator in Kirshner, H., et al (2013) J. Microsc. 
%   (https://bigwww.epfl.ch/algorithms/psfgenerator/). A measured
%   instrument psf may be a better option -YMMV.
%   For the Gaussian function, sigma (px) is calulated from FWHM (micron) 
%   set in section 1. The kernel size is calculated and dimensions are odd
%   so that the Gaussian function is centred on a pixel.
psfFwhm=[fwhmXY fwhmXY fwhmZ];      
psfSigma=[psfFwhm(1:2)/dxy psfFwhm(3)/dz]/2.355;    
psfSz=round(psfSigma*6);                        
isEven=mod(psfSz,2)==0;                         
psfSz(isEven)=psfSz(isEven)+1;
psf=fspecial3("gaussian",psfSz,psfSigma);

%   4b. Deconvolve data and inspect result: Here, a Richardson-Lucy 
%   deconvolution algorithm is used to improve signal-to-noise and data is 
%   padded to limit edge effects. Some ringing due to incorrect PSF data
%   is not so serious as these data are only used for skeleton
%   construction...
iterations=7;                                         
padSz=psfSz;      
padData=padarray(data,padSz,"both","replicate");
deconData=deconvlucy(padData,psf,iterations);
deconData=deconData(padSz(1)+1:end-padSz(1),padSz(2)+1:end-padSz(2),...
                    padSz(3)+1:end-padSz(3)); %   strip padding

%   4c. Create data-dependent masks for the bath, cell interior and t-tubules: 
%   Inspect output carefully to ensure that the masks reflect the input data. 
%   Here, the initial ttMask obtained from adaptive thresholding and threshold 
%   values are data sensitive -adjust as needed prior to skeletonisation. 
scaledData=imgaussfilt3(data,1./dxy);
scaledData=rescale(scaledData);
cellMask=scaledData<graythresh(scaledData);
bathMask=~cellMask;
cellMask=scaledData<graythresh(scaledData(cellMask(:)));
cellMask=imerode(cellMask,strel("disk",round(1/dxy)));
bathMask=scaledData>graythresh(scaledData(bathMask(:)));
bathMask=imerode(bathMask,strel("disk",round(1/dxy)));
ttMask=rescale(deconData);
nSz=round([1.8/dxy 1.8/dxy 1.0/dz]) ; nSz(mod(nSz,2)==0)=nSz(mod(nSz,2)==0)+1;
localT=adaptthresh(ttMask,0.4,"NeighborhoodSize",nSz);
ttMask=ttMask>localT;
unitSphere=strel("sphere",1);
ttMask=imclose(ttMask,unitSphere);              % (optional) tidy mask, 
ttMask=imerode(ttMask,unitSphere);              % e.g. close, erode, clean
ttMask=bwmorph3(ttMask,"clean");                
ttMask=and(ttMask,cellMask);                       % mask off bath volume
ttSkel=bwskel(padarray(ttMask,[1 1 1],"both"));    % make 3d skeleton
ttSkel=ttSkel(2:end-1,2:end-1,2:end-1);
ttSkel=bwmorph3(ttSkel,"clean");  


%%  5.  Analyze data
%   5a. Asymmetric blurring of data to remove effect of asymmetric PSF
%   (aka "spherising" data). Normalise data to bath signal.
spherSigma=[psfFwhm(3)/dz psfFwhm(3)/dz psfFwhm(1)/dxy]/2.355;
spherData=imgaussfilt3(data.*cellMask,spherSigma);
bathData=imgaussfilt3(data.*bathMask,spherSigma);
spherData=spherData/median(bathData(bathMask));
ttData=spherData.*ttSkel;
fprintf("Mean t-tubule radius = %1.4f (in weighted R units)\n",mean(ttData(ttSkel)));
figure(1), hold off, histogram(ttData(ttSkel),0:0.01:0.4);

%  5b.Correct the skeleton for curving and branching (called W in SC 1999). 
%  To convert to physical values you need to solve equation 4 given in 
%  Soeller C, Cannell MB. (1999) Circulation research. 84:266-275 
%  for the actual PSF used. 
ttm=zeros([9 9 9]) ; ttm(5,:,5)=1 ; % make tt model
Sphr=getnhood(strel("sphere",2)).*1.0;
norm=convn(ttm,Sphr,"same");   % probe tt model with Sphr to get W for straight tt
norm=max(norm(:)).*(fwhmXY/fwhmZ);
wSkel=convn(ttSkel,Sphr,"same")./norm; % probe skeleton geoemtry with sphere
corData=ttData./wSkel;
figure(1), hold on, histogram(corData(ttSkel),0:0.01:0.4);
legend('wR','R')
fprintf("Mean t-tubule radius = %1.4f (in R/Rinf, with branch weight correction)\n",mean(corData(ttSkel)));
vol=volshow(wSkel,Colormap=colmap); % show skeleton weights
