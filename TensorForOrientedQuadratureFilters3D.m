function T = TensorForOrientedQuadratureFilters3D(im,orient)
%%  Tensor - calculate 3D tensor of the image
%   
%   REFERENCE:
%       Hans Knutsson, 
%       Representing local structure using tensors, 
%       The 6th Scandinavian Conference on Image Analysis, Oulu, Finland, 
%       244–251, June 1989
%
%   INPUT:
%       im - 3D gray image
%       s  - Gaussian kernel sigma
%
%   OUTPUT:
%       T  - 3D tensor
%
%   HELP:
%       T = SUMi(|qi|(nini' - alphaI) 
%       where: 
%           qi - response of i-th quadrature filter in the spatial domain
%           alpha = 1/(m-1), m - dimensionality of T (2D -> m=4; 3D -> m=6)
%           I - identity tensor
%
%   USAGE:
%      T = TensorForOrientedQuadratureFilters3D(im,5);
%
%   AUTHOR:
%       Boguslaw Obara, http://boguslawobara.net/
%
%   VERSION:
%       0.1 - 16/07/2012 First implementation

%% Tensor
I = eye(3,3);
T = zeros(size(im,1),size(im,2),size(im,3),9);
for i=1:size(orient,1)
   nnT = kron(orient(i,:),orient(i,:)');
   %nnTI = (nnT-1/(4-1).*I); %alpha = 1/(m-1) where m = size of T
   %nnTI = (nnT-(1/(2-1)).*I); %alpha = 1/(m-1) where m = dimensionality of T    
   %nnTI = (nnT-(1/(length(ang)-1)).*I); %alpha = 1/(m-1) where m = nr of filters
   %nnTI = (nnT-(1/(3-1)).*I); %alpha = 1/(m+1) where m = dimensionality of T (A FRAMEWORK FOR ESTIMATION OF ORIENTATION AND VELOCITY)     
   
   nnTI = nnT;
   T(:,:,:,1) = T(:,:,:,1) + im(:,:,:,i)*nnTI(1,1);
   T(:,:,:,2) = T(:,:,:,2) + im(:,:,:,i)*nnTI(1,2);
   T(:,:,:,3) = T(:,:,:,3) + im(:,:,:,i)*nnTI(1,3);
   
   T(:,:,:,4) = T(:,:,:,4) + im(:,:,:,i)*nnTI(2,1);
   T(:,:,:,5) = T(:,:,:,5) + im(:,:,:,i)*nnTI(2,2);
   T(:,:,:,6) = T(:,:,:,6) + im(:,:,:,i)*nnTI(2,3);
   
   T(:,:,:,7) = T(:,:,:,7) + im(:,:,:,i)*nnTI(3,1);
   T(:,:,:,8) = T(:,:,:,8) + im(:,:,:,i)*nnTI(3,2);
   T(:,:,:,9) = T(:,:,:,9) + im(:,:,:,i)*nnTI(3,3);    
end
%% 
end