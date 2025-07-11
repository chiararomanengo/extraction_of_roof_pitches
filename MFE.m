function mfe=MFE(xyz,dist)

% Input: a point cloud xyz and a vector containing the Euclidean distances
% Output: the mean fitting error, that is a percentage error computed as
% the mean of distances in dist weighted with respect to the diagonal of 
% the axis aligned bounding box of xyz

% Computation of the diagonal of the axis-aligned bounding box
base=max(xyz(:,1))-min(xyz(:,1));
h=max(xyz(:,2))-min(xyz(:,2));
diag=sqrt(base^2+h^2);
h=max(xyz(:,3))-min(xyz(:,3));
Fin=sqrt(diag^2+h^2);

% Mean fitting error
mfe=mean(dist)/Fin;

end