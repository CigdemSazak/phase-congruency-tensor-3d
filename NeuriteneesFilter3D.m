function [imf,L1,L2,L3,Lmin] = NeuriteneesFilter3D(im,sigma,saveFlag,name)
%%  NeuriteneesFilter - computing the neuritenees
%   
%   REFERENCE:
%       E. Meijering et al. Design and Validation of a Tool for Neurite 
%       Tracing and Analysis in Fluorescence Microscopy Images
%       Cytometry Part A, 58A, 167–176, 2004
%
%   INPUT:
%       im      - gray image
%       sigma   - sigma factor
%
%   OUTPUT:
%       imf     - filtered image
%
%   USAGE:
%
%   AUTHOR:
%       Boguslaw Obara, http://boguslawobara.net/
%
%   VERSION:
%       0.1 - 25/04/2010 First implementation

%% TENSOR

T = TensorForOrientedQuadratureFilters3D(im,sigma);

%% Eigen Matrix - values and vectors

[L1,L2,L3,V1,V2,V3,V4,V5,V6,V7,V8,V9] ...
    = EigenMatrix3x3M(T(:,:,:,1),T(:,:,:,2),T(:,:,:,3),...
                        T(:,:,:,4),T(:,:,:,5),T(:,:,:,6),...
                        T(:,:,:,7),T(:,:,:,8),T(:,:,:,9));

%% Modified Hessian
alfa = -1/3;

L1p = L1 + alfa/2.*L2 + alfa/2.*L3;
L2p = L2 + alfa/2.*L1 + alfa/2.*L3;
L3p = L3 + alfa/2.*L1 + alfa/2.*L2;
L1 = L1p;
L2 = L2p;
L3 = L3p;

%% Sort L1s > L2s> L3s
tic;
[L1,L2,L3,V1,V2,V3,V4,V5,V6,V7,V8,V9] = EigenSort3x3M(L1,L2,L3,V1,V2,V3,V4,V5,V6,V7,V8,V9);
toc;
%% Neuriteness
Lmin = min(L1(:));
imf = zeros(size(L1));
imf(L1<0) = L1(L1<0)./Lmin;

%% Save
if(saveFlag == 1);
imf = Normalize(imf);
s2D = strcat('PCTNeuriteNew', name, '.png');
s3D = strcat('PCTNeurite', name, '.tif');
saveFunction(imf,'2D',s2D);
% saveFunction(imf,'3D',s3D);
end
end
