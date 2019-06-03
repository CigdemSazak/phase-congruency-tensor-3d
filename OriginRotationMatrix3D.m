function [ux,uy,uz,R] = OriginRotationMatrix3D(ux,vx) 
%%
if nargin<2
    vx = [1 0 0];
    %vy = [0 1 0];
    %vz = [0 0 1];
end
%%
vx = vx/norm(vx);
%vy = vy/norm(vy);
%vz = vz/norm(vz);
%ux = [1 1 1];
%uz = cross(vx,ux);
%uy = cross(uz,ux);
uy = cross(ux,vx);
uz = cross(ux,uy);
ux = ux/norm(ux);
uy = uy/norm(uy);
uz = uz/norm(uz);
%% Rotation Matrix
if abs(ux) == abs(vx)
    R = eye(3,3);
    R = sign(ux(1))*R;
else
    R = [ux; uy; uz];
end