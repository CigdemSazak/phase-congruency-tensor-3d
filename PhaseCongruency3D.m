function [PC,EO,T,pcSum,orients] = PhaseCongruency3D(im,nscale,norient,minWaveLength,mult,sigmaOnf,k,cutOff,g,noiseMethod)
%% Args
% Get arguments and/or default values    
% [im, nscale, norient, minWaveLength, mult, sigmaOnf, ...
%                    k, cutOff, g, noiseMethod] = checkargs(varargin(:));
%    pc         - Cell array of phase congruency images (values between 0 and 1)   
%                 for each orientation
%%
epsilon         = .0001;          % Used to prevent division by zero.
[rows,cols,planes] = size(im);
imagefft = fftn(im);              % Fourier transform of image

[p,ps,dm] = SurfacesSpiralPoints3D(norient);% norient = 7, 40
% ps sphere spiral points
orients = ps;
norient = size(orients,1);

zero = zeros(rows,cols,planes);
EO = cell(nscale, norient);       % Array of convolution results.  
PC = cell(norient,1);
% covx2 = zero;                     % Matrices for covariance data
% covy2 = zero;
% covxy = zero;

EnergyV = zeros(rows,cols,planes,3);     % Matrix for accumulating total energy
                                  % vector, used for feature orientation
                                  % and type calculation

pcSum = zeros(rows,cols,planes);                                      
% Pre-compute some stuff to speed up filter construction

% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols; 
end

if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows; 
end

if mod(planes,2)
    zrange = [-(planes-1)/2:(planes-1)/2]/(planes-1);
else
    zrange = [-planes/2:(planes/2-1)]/planes;	
end

[x,y,z] = meshgrid(xrange, yrange, zrange);

[theta,phi,radius] = cart2sph(x,y,z);
% radius = sqrt(x.^2 + y.^2 + z.^2);% Matrix values contain *normalised* radius from centre.
% radius2 = sqrt(x.^2 + y.^2); 
% theta = atan2(-y,x);              % Matrix values contain polar angle.
%                                   % (note -ve y is used to give +ve
%                                   % anti-clockwise angles)
% phi = atan2(radius2,z);

radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.
phi    = ifftshift(phi);
radius(1,1,1) = 1;                % Get rid of the 0 radius value at the 0
                                  % frequency point (now at top-left corner)
                                  % so that taking the log of the radius will 
                                  % not cause trouble.

% THIS IS POLAR TO CARTESIAN TRANSFORMATION                                  
%x = sin(theta).*cos(phi);   % x
%y = sin(theta).*sin(theta); % y
%z = cos(theta);             % z
[x,y,z] = sph2cart(theta,phi,radius);

%clear x; clear y; 
clear theta;    % save a little memory

% Filters are constructed in terms of two components.
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.

% Construct the radial filter components...
% First construct a low-pass filter that is as large as possible, yet falls
% away to zero at the boundaries.  All log Gabor filters are multiplied by
% this to ensure no extra frequencies at the 'corners' of the FFT are
% incorporated as this seems to upset the normalisation process when
% calculating phase congrunecy.
lp = lowpassfilter3D([rows,cols,planes],.45,15);   % Radius .45, 'sharpness' 15

logGabor = cell(1,nscale);

for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
                                          % back to zero (undo the radius fudge).
    %figure; imagesc(fftshift(logGabor{s}));
end
imx = zeros(size(im));
for o = 1:norient                    % For each orientation...
    % Construct the angular filter spread function
    %[thetao,phio,radiuso] = cart2sph(orients(o,1),orients(o,1),orients(o,2));

    
    ux = orients(o,:); 
    %ux = [1 1 1];
    vx = [1 0 0]; %vy = [0 1 0]; vz = [0 0 1];
    [ux,uy,uz,R] = OriginRotationMatrix3D(ux,vx); 
    %[ a*r11 + b*r21 + c*r31, a*r12 + b*r22 + c*r32, a*r13 + b*r23 + c*r33]
    xp = x*R(1,1) + y*R(1,2) + z*R(1,3);
    yp = x*R(2,1) + y*R(2,2) + z*R(2,3);
    zp = x*R(3,1) + y*R(3,2) + z*R(3,3);
    xo = ux(1)*R(1,1) + ux(2)*R(1,2) + ux(3)*R(1,3);
    yo = ux(1)*R(2,1) + ux(2)*R(2,2) + ux(3)*R(2,3);
    zo = ux(1)*R(3,1) + ux(2)*R(3,2) + ux(3)*R(3,3);        

 
    % For each point in the filter matrix calculate the angular distance from
    % the specified filter orientation.  To overcome the angular wrap-around
    % problem sine difference and cosine difference values are first computed
    % and then the atan2 function is used to determine angular distance.
    % THIS IS ROTATION !!!

    %ds = sintheta * cos(angl1) - costheta * sin(angl1);    % Difference in sine.
    %dc = costheta * cos(angl1) + sintheta * sin(angl1);    % Difference in cosine.
    %dtheta = abs(atan2(ds,dc));                            % Absolute angular distance.

    xo = ones(size(im))*xo; 
    yo = ones(size(im))*yo; 
    zo = ones(size(im))*zo;
    dtheta = VectorVectorAngle3D(xp,yp,zp,xo,yo,zo);
    %dtheta = BOVectorVectorAngle3D([xp yp zp],[xo yo zo]);
    %dtheta = abs(dtheta);
    
    % Scale theta so that cosine spread function has the right wavelength and clamp to pi    
    %dtheta = min(dtheta*norient/2,pi);
    % Ratio = (1/3)*sqrt(3)*a-1/2*a = ((2*sqrt(3)-3)/6)*a;
    % A1 = pi-3/2*sqrt(3); A2 = 3*pi; Ratio = A1/A2;
    Ratio = power(((2*sqrt(3)-3)/6)*1,1/3);
    norient2 = norient-norient*Ratio;
    %norient2 = norient;
    dtheta = min(dtheta*norient2/2,pi); 
    % The spread function is cos(dtheta) between -pi and pi.  We add 1,
    % and then divide by 2 so that the value ranges 0-1
    spread = (cos(dtheta)+1)/2;  
    
    % figure; 
    % [xg,yg,zg] = meshgrid(1:size(im,2),1:size(im,1),1:size(im,3));
    % p = patch(isosurface(xg,yg,zg,ifftshift(spread)>0,0.5));
    % set(p,'FaceColor','red','EdgeColor','none');
    % daspect([1 1 1]); view(3); 
    % box on; grid on;
    % camlight; lighting gouraud;
    % xlim([1 size(im,2)]) 
    % ylim([1 size(im,1)]) 
    % zlim([1 size(im,3)]) 
    % title([num2str(ux(1)) '; ' num2str(ux(2)) '; ' num2str(ux(3))])
    imx(ifftshift(spread)>0) = 1;

    sumE_ThisOrient   = zero;          % Initialize accumulator matrices.
    sumO_ThisOrient   = zero;       
    sumAn_ThisOrient  = zero;      
    Energy            = zero;      

    for s = 1:nscale,                  % For each scale...
        filter = logGabor{s} .* spread;      % Multiply radial and angular
                                             % components to get the filter. 

        % Convolve image with even and odd filters returning the result in EO
        EO{s,o} = ifftn(imagefft .* filter);      

        An = abs(EO{s,o});                         % Amplitude of even & odd filter response.
        sumAn_ThisOrient = sumAn_ThisOrient + An;  % Sum of amplitude responses.
        sumE_ThisOrient = sumE_ThisOrient + real(EO{s,o}); % Sum of even filter convolution results.
        sumO_ThisOrient = sumO_ThisOrient + imag(EO{s,o}); % Sum of odd filter convolution results.
        %% ATTENTION !!!!
        % the orientation is calculated for edges and not for linear
        % objects, you can use:
        %sumO_ThisOrient = sumO_ThisOrient + real(EO{s,o}); % Sum of even filter convolution results.

        % At the smallest scale estimate noise characteristics from the
        % distribution of the filter amplitude responses stored in sumAn. 
        % tau is the Rayleigh parameter that is used to describe the
        % distribution.
        if s == 1 
            if noiseMethod == -1     % Use median to estimate noise statistics
                tau = median(sumAn_ThisOrient(:))/sqrt(log(4));   
            elseif noiseMethod == -2 % Use mode to estimate noise statistics
                tau = rayleighmode(sumAn_ThisOrient(:));
            end
            maxAn = An;
        else
            % Record maximum amplitude of components across scales.  This is needed
            % to determine the frequency spread weighting.
            maxAn = max(maxAn,An);   
        end
    end                                       % ... and process the next scale
    
    % Accumulate total 3D energy vector data, this will be used to
    % determine overall feature orientation and feature phase/type
    EnergyV(:,:,:,1) = EnergyV(:,:,:,1) + sumE_ThisOrient;
    %EnergyV(:,:,:,2) = EnergyV(:,:,:,2) + cos(angl)*sumO_ThisOrient;
    %EnergyV(:,:,:,3) = EnergyV(:,:,:,3) + sin(angl)*sumO_ThisOrient;

    % Get weighted mean filter response vector, this gives the weighted mean
    % phase angle.
    XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;   
    MeanE = sumE_ThisOrient ./ XEnergy; 
    MeanO = sumO_ThisOrient ./ XEnergy; 

    % Now calculate An(cos(phase_deviation) - | sin(phase_deviation)) | by
    % using dot and cross products between the weighted mean filter response
    % vector and the individual filter response vectors at each scale.  This
    % quantity is phase congruency multiplied by An, which we call energy.

    for s = 1:nscale,       
        E = real(EO{s,o}); O = imag(EO{s,o});    % Extract even and odd
                                                 % convolution results.
        Energy = Energy + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
    end

    %% Automatically determine noise threshold
    %
    % Assuming the noise is Gaussian the response of the filters to noise will
    % form Rayleigh distribution.  We use the filter responses at the smallest
    % scale as a guide to the underlying noise level because the smallest scale
    % filters spend most of their time responding to noise, and only
    % occasionally responding to features. Either the median, or the mode, of
    % the distribution of filter responses can be used as a robust statistic to
    % estimate the distribution mean and standard deviation as these are related
    % to the median or mode by fixed constants.  The response of the larger
    % scale filters to noise can then be estimated from the smallest scale
    % filter response according to their relative bandwidths.
    %
    % This code assumes that the expected reponse to noise on the phase congruency
    % calculation is simply the sum of the expected noise responses of each of
    % the filters.  This is a simplistic overestimate, however these two
    % quantities should be related by some constant that will depend on the
    % filter bank being used.  Appropriate tuning of the parameter 'k' will
    % allow you to produce the desired output. 

    if noiseMethod >= 0     % We are using a fixed noise threshold
        T = noiseMethod;    % use supplied noiseMethod value as the threshold
    else
        % Estimate the effect of noise on the sum of the filter responses as
        % the sum of estimated individual responses (this is a simplistic
        % overestimate). As the estimated noise response at succesive scales
        % is scaled inversely proportional to bandwidth we have a simple
        % geometric sum.
        totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult));

        % Calculate mean and std dev from tau using fixed relationship
        % between these parameters and tau. See
        % http://mathworld.wolfram.com/RayleighDistribution.html
        EstNoiseEnergyMean = totalTau*sqrt(pi/2);        % Expected mean and std
        EstNoiseEnergySigma = totalTau*sqrt((4-pi)/2);   % values of noise energy

        T =  EstNoiseEnergyMean + k*EstNoiseEnergySigma; % Noise threshold
    end

    % Apply noise threshold,  this is effectively wavelet denoising via
    % soft thresholding.
    Energy = max(Energy - T, 0);         

    % Form weighting that penalizes frequency distributions that are
    % particularly narrow.  Calculate fractional 'width' of the frequencies
    % present by taking the sum of the filter response amplitudes and dividing
    % by the maximum amplitude at each point on the image.   If
    % there is only one non-zero component width takes on a value of 0, if
    % all components are equal width is 1.
    width = (sumAn_ThisOrient./(maxAn + epsilon) - 1) / (nscale-1);    

    % Now calculate the sigmoidal weighting function for this orientation.
    weight = 1.0 ./ (1 + exp( (cutOff - width)*g)); 

    % Apply weighting to energy and then calculate phase congruency
    PC{o} = weight.*Energy./sumAn_ThisOrient;   % Phase congruency for this orientation

    pcSum = pcSum + PC{o};        
end  % For each orientation

end
%% End