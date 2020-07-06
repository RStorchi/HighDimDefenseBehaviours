function[Xest] = ls_triangulate(x,P)
%this is the single point LS linear triangulation for an arbitrary number
%of cameras
%INPUT: 
%x: homogeneous 2D coordinates
%P: camera matrices (array of cells, one cell per camera matrix)
%OUTPUT: 
%Xest: 3D estimate in normalized homogeneous coordinates 
 
Ncam = length(P);
A = zeros(2*Ncam,4);
for n = 1:Ncam
    A(2*(n-1)+1,:) = x(1,n)*P{n}(3,:)-P{n}(1,:);
    A(2*n,:) = x(2,n)*P{n}(3,:)-P{n}(2,:);
end
[U,S,V] = svd(A); Xest = V(:,end);
Xest = Xest/Xest(4);