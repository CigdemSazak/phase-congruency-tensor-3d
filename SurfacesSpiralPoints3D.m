function [p,ps,dm] = SurfacesSpiralPoints3D(n)
%% Origin
pn = [1 0 0];
pn = pn/norm(pn);
%% Spiral points
p = sphere_spiralpoints(1,zeros(3,1),n)';
%% Angle & Distance
a1 = zeros(size(p,1),size(p,1));
a2 = zeros(size(p,1),1);
d = zeros(size(p,1),size(p,1));
for i=1:size(p,1)
    for j=1:size(p,1)
        if i~=j
            a1(i,j) = VectorVectorAngle3D(p(i,1),p(i,2),p(i,3),p(j,1),p(j,2),p(j,3));
            d(i,j) = norm(p(i,:)-p(j,:));
        end
    end
    a2(i) = VectorVectorAngle3D(p(i,1),p(i,2),p(i,3),pn(1),pn(2),pn(3));
end
%% Min Distance
dm = [];
for i=1:size(p,1)
    idx = 1:size(p,1);
    idx(i) = [];
    dm(i) = min(d(i,idx));
end
%% Prune half of the points
idx = find(a2>pi/2);
%% Prune 2pi
[xi,yi] = find(a1==pi);
if ~isempty(xi)
    xi(ismember(xi,idx)) = [];
    xi = xi(1:length(xi)/2);
    idx = [idx; xi];
end
%% Prune
ps = p;
ps(idx,:) = [];
dm(idx) = [];
%% End
end