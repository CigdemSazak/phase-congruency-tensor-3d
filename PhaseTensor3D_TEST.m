%% Clear
clc; clear all;
%close all;
%% synthetic Data
im = syntheticData3D(2);
%% Normalize
im = Normalize(im);
%% Phase
nscale          = 3;     % Number of wavelet scales.    
norient         = 10;    % Number of filter orientations.
minWaveLength   = 3;     % Wavelength of smallest scale filter.    
mult            = 2.1;   % Scaling factor between successive filters.    
sigmaOnf        = 0.45;  % Ratio of the standard deviation of the
                         % Gaussian describing the log Gabor filter's
                         % transfer function in the frequency domain
                         % to the filter center frequency.    
k               = 5;     % No of standard deviations of the noise
                         % energy beyond the mean at which we set the
                         % noise threshold point. 
cutOff          = 0.5;   % The fractional measure of frequency spread
                         % below which phase congruency values get penalized.
g               = 10;    % Controls the sharpness of the transition in
                         % the sigmoid function used to weight phase
                         % congruency for frequency spread.                      
noiseMethod     = -1;    % Choice of noise compensation method. 
%% Phase Congruency 3D
[PC,EO,T,pcSum,orients] = PhaseCongruency3D(im,nscale,norient,minWaveLength,mult,sigmaOnf,k,cutOff,g,noiseMethod,saveFlag,name);
%%
pcSum = Normalize(pcSum);
imq = zeros(size(im,1),size(im,2),size(im,3),size(PC,1));
for i=1:size(PC,1)
    imq(:,:,:,i) = PC{i,1};
end
imq = Normalize(imq);
%% Phase Vesselness 3D
alpha = 0.5; beta = 0.5; c = 15;
Vmax = PhaseVesselness3D(imq,orients,alpha,beta,c,saveFlag,name);

%% Phase Neuriteness 3D
 sigma = [2 2 2] ; 
[imf,L1,L2,L3,Lmin] = NeuriteneesFilter3D(imq,sigma,saveFlag,name);
%% Display
figure, imagesc(max(im,[],3)); colormap jet; axis equal; axis tight; axis off;
figure, imagesc(max(imf,[],3)); colormap jet; axis equal; axis tight; axis off;
figure, imagesc(max(Vmax,[],3));colormap jet; axis equal; axis tight; axis off;
