function [coord, maxCoord]=PlanesHT0_search(X)

% Input: a point cloud X
% Output: parameter of the best fitting plane (coord) and the maximum value
% of the accumulator matrix (maxCoord)

% Inizialization of the parameter space of the plane in the Hesse form 


rho= max(sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2));
rho_min= min(sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2));
c=(pi*2)/180;
step=((rho+rho_min)/2)*sin(2*c);

a=0:step:rho;
b=0:3:357;
b=(pi*b)/180;
c=0:2:178;
c=(pi*c)/180;

soglia=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
soglia(find(soglia==0))=10^(-2);

Na=numel(a);
Nb=numel(b);
Nc=numel(c);

% Inizialization of the accumulator matrix

H=zeros(Na,Nb,Nc);

% Estimation of the accumulator matrix

for j=1:Na
    for k=1:Nb
        for h=1:Nc
            H(j,k,h)=numel(find(abs(a(j)-X(:,1)*cos(b(k))*sin(c(h))-X(:,2)*sin(c(h))*sin(b(k))-X(:,3)*cos(c(h)))<soglia/500));
        end
    end
end

maxCoord=max(max(max(H)));

% Searching for the coordinates corresponding the maximum value of the
% accumulator

coord=[];
for j=1:Na
    for k=1:Nb
        for h=1:Nc
            if  H(j,k,h)==maxCoord
                coord=[coord; a(j) b(k) c(h)];
            end            
        end
    end
end


end