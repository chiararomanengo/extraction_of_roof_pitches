function [output, mfe,thetaf, phif] = SearchPlane(xyz)

% Input: a point cloud xyz
% Output: a vector output containing the parameter a, b, c, d of the plane
% in the form ax+by+cz+d=0
% 
try
    
    % Translation of the point cloud: barycenter in the origin
    
    trX=mean(xyz(:,1));
    trY=mean(xyz(:,2));
    trZ=mean(xyz(:,3));
    
    xyz(:,1)=xyz(:,1)-trX;
    xyz(:,2)=xyz(:,2)-trY;
    xyz(:,3)=xyz(:,3)-trZ;
    
    % Application of the Hough Transform to find the best fitting plane in
    % the Hesse form: coord contains the parameter rho, theta and phi of 
    % the equation cos(theta)sin(phi)x + sin(phi)sin(theta)y + cos(phi) +
    % - rho=0
    
    [coord, maxCoord]=PlanesHT0_search(xyz);
    
    maxCoord
    % In case, the number of fitting plane is higher than 1, the following
    % step fixes the one with lower approximation error (MFE)
    
    mAus=1;
    for j=1:size(coord,1)
        rho=coord(j,1);
        theta=coord(j,2);
        phi=coord(j,3);
        if abs(cos(phi))>10^(-3)
            x_piano = linspace(min(xyz(:,1)),max(xyz(:,1)),200);
            y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),200);
            piano = zeros(length(x_piano)*length(y_piano),3);
            for i=1:length(x_piano)
                for k=1:length(y_piano)
                    piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                    piano( (i-1)*length(x_piano)+k,2) = y_piano(k);
                    piano( (i-1)*length(x_piano)+k,3) = (rho-x_piano(i)*cos(theta)*sin(phi)-y_piano(k)*sin(phi)*sin(theta))/cos(phi);
                end
            end
        else
            if abs(cos(theta))>10^(-3)
                y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),200);
                z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),200);
                piano = zeros(length(y_piano)*length(z_piano),3);
                for i=1:length(y_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(y_piano)+k,1) = (rho-y_piano(i)*sin(theta)*sin(phi))/(sin(phi)*cos(theta));
                        piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                        piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                    end
                end
            else
                x_piano=linspace(min(xyz(:,1)),max(xyz(:,1)),200);
                z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),200);
                piano = zeros(length(x_piano)*length(z_piano),3);
                for i=1:length(x_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                        piano( (i-1)*length(x_piano)+k,2) = (rho)/(sin(phi)*sin(theta));
                        piano( (i-1)*length(x_piano)+k,3) = z_piano(k);
                    end
                end
            end
        end
        
        [I,dist] = knnsearch(piano,xyz);
        m=MFE(xyz,dist);
        
        if m<mAus
            mAus=m;
            mfe=m;
            output=[cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi) cos(theta)*sin(phi)*trX+sin(phi)*sin(theta)*trY+cos(phi)*trZ-rho]; %cos(theta)*sin(phi)*trX+sin(phi)*sin(theta)*trY+cos(phi)*trZ
            RHO=rho;
            thetaf=theta;
            phif=phi;
        end
    end
    

    
catch
    mfe=NaN;
    output=[];
%     mBB=[];
%     roof=[];
    thetaf=[];
    phif=[];
    return
end



end