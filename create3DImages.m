
%%  create3DImages - 3D Synthetic image
%   
%   AUTHOR:
%       Cigdem Sazak,
%
%   VERSION:
%       0.1 - 16/07/2016 First implementation
function im = create3DImagesingle
%% Set-Up
angle = mod(90,180);
vesselRadius = 4;%gives thickness of 9
theta = angle/2;
sep =1;
blobRadius = 4;
im = zeros(100,100,100);

%% Incoming Vessel
im(50,50,1:50) = 1;
im(50,1:99,50) = 1;
im(1:99,50,50) = 1;
%% Outgoing Vessels
dx = linspace(0,50,5100);
%Upper
dy = dx*tand(theta);
x = dx+50;
y = dy+50; 
z = dx+50;x(y>100)=[]; z(y>100)=[]; y(y>100)=[];
index = sub2ind(size(im),round(y),round(x),round(z));%indices of the line
im(index) = 1;
% Lower
dy = dx*tand(-theta);
x = dx+50;
y = dy+50; 
z = dx+50; x(y<0.5)=[]; z(y<0.5)=[]; y(y<0.5)=[]; 
index = sub2ind(size(im),round(y),round(x),round(z));%indices of the line
im(index) = 1;

%% create blob
gtB = zeros(100,100,100);
gtB(50,round(45-(sep*max(vesselRadius,blobRadius))),25) = 1;
gtB = bwdist(gtB)<blobRadius;
im = max(im,gtB);
%%
im = bwdist(im)<4;
im= double(im);
im = imnoise(im,'salt & pepper');
im = imgaussfilt3(im,2);
% figure,imshow3D(im);
% %%
% % figure,imagesc(squeeze(max(im,[],1)));colormap gray;
figure; isosurface(im);
axis equal;
% axis([1 100 1 100 1 100]);
% view(30,45);
grid on
camlight
colormap copper
xlabel('X');
ylabel('Y');
zlabel('Z');
% print('Synthetic3D','-dpng')
%%
% gt = zeros(102,100,100);
%  im =imread('3dSynthetic.png');
%  im = double(im);
%  im = imbinarize(rgb2gray(im));
%  B = bwdist(im);
% %  B = bwmorph(~im, 'skel');
%  gt(:,:,50) = B;
%  im = bwdist(gt);
% % figure, histogram(gt);
% %  figure, imshow3D(gt);
% % im = imcomplement(im);
%  figure,imagesc(max(im,[],3));colormap gray;
%  
%  figure, 
%  isosurface(im); axis equal;
%  axis([1 100 1 100 1 100]);
% % view(30,45);
% grid on
% camlight
% colormap copper
end
