clc; clear all; close all

%% Converting .las into .txt -- > x y z nx ny nz

dirname=strcat('../vlp/');
Files_cell=dir(dirname);
Files_cell(1:2)=[];
mfeFINALE=[];

%%
for i=1:size(Files_cell,1)
    
    FoldNames=Files_cell(i).name;
    FoldNames=strcat(dirname, FoldNames);
    Files=dir(FoldNames);
    
    if size(Files,1)>2
        Files(1:2)=[];
        FileNames=Files(1).name;
        
        aus=double(string(Files_cell(i).name));
        INDEX_ID=aus;
        
        FileNames=strcat(FoldNames, '/', FileNames);
        lasReader = lasFileReader(FileNames);
        [ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","Classification");
        xyz=ptCloud.Location;
        labels = label2rgb(pointAttributes.Classification);
        colorData = reshape(labels,[],3);
        
        if size(xyz,1)>50           
            
             dirName=strcat('buildings_Segm/building',  int2str(aus));
             [status, msg, msgID] = mkdir(dirName);
            nameFile=strcat(dirName, '/building',  int2str(aus), '.las');
            
            
            fileID = fopen(nameFile,'w');
            fprintf(fileID,'%14.12f %14.12f %14.12f \n', xyz(1:(size(xyz,1)-1),:)');
            fprintf(fileID,'%14.12f %14.12f %14.12f \n', xyz(size(xyz,1),:)');
            fclose(fileID);
            
            outFile=strcat(dirName, '/building_norm',  int2str(aus), '.xyz');
            command = strcat('cgal-scaline-normals.exe -i',{' '}, FileNames,' -o ',{' '}, outFile);
            command=command{1,1};
            status = system( sprintf(command) );
            
            outFile=strcat(dirName, '/building_norm',  int2str(aus), '_angle_and_flag.xyz');
%             outFile=strcat(dirName, '/building_norm',  int2str(aus), '_nothing.xyz');
            [x, y, z]=xyzread(outFile);
            nx=x(2:2:end);
            x=x(1:2:end);
            ny=y(2:2:end);
            y=y(1:2:end);
            nz=z(2:2:end);
            z=z(1:2:end);
            
            xyz=[x y z];
            normals=[nx ny nz];
            
            trX=mean(xyz(:,1));
            trY=mean(xyz(:,2));
            trZ=mean(xyz(:,3)); % ricorda traslazione!!!
            xyz(:,1)=xyz(:,1)-trX;
            xyz(:,2)=xyz(:,2)-trY;
            xyz(:,3)=xyz(:,3)-trZ;
            
            xyz_norm=[xyz normals];
            nameFile=strcat(dirName, '/building',  int2str(aus), '.txt');
            fileID = fopen(nameFile,'w');
            fprintf(fileID,'%14.12f %14.12f %14.12f %14.12f %14.12f %14.12f \n', xyz_norm(1:(size(xyz_norm,1)-1),:)');
            fprintf(fileID,'%14.12f %14.12f %14.12f %14.12f %14.12f %14.12f \n', xyz_norm(size(xyz_norm,1),:)');
            fclose(fileID);
            
            nameFile=strcat(dirName, '/building',  int2str(aus), '.txt');
            outFile=strcat(dirName, '/building_output',  int2str(aus), '.ply');
            num=size(xyz,1)*0.01/4;
            
            if num<10
                num=10;
            end
             
            command = strcat('fit_region_growing-fit.exe -i',{' '}, nameFile,' -o ',{' '}, outFile, ' -a 10 -s',{' '}, int2str(num));
            command=command{1,1};
            status = system( sprintf(command) );
                      
            obj=pcread(outFile);
            
            xyz=obj.Location;
            color=obj.Color;
            
            mBB_max_X=max(xyz(:,1));
            mBB_min_X=min(xyz(:,1));
            mBB_max_Y=max(xyz(:,2));
            mBB_min_Y=min(xyz(:,2));
            
            control=double(color);
%             control=unique(control,'rows');
            
            segment={};
            aus=[0 0 0];
            for s=1:size(control,1)
                if norm(aus-control(s,:))>0
                    segment=[segment; xyz(find(sum(color==control(s,:),2)>2),:)];
                    aus=double(control(s,:));
                end
            end
            
            Final=zeros(size(segment,1),4);
            angles=zeros(size(segment,1),2);
            lab_vert=[];
            angles_orig=zeros(size(segment,1),2);
            for s=1:size(segment,1)
                segm_xyz=segment{s,1};
                [output, mfe, theta, phi]=SearchPlane(segm_xyz);
                mfeFINALE=[mfeFINALE; mfe];
                Final(s,:)=output;

                angles(s,1)=theta*180/pi;
                angles(s,2)=phi*180/pi;
                angles_orig(s,1)=theta*180/pi;
                angles_orig(s,2)=phi*180/pi;

                if angles(s,2)>90
                    angles(s,2)=180-angles(s,2);
                    if  angles(s,1)<180
                        angles(s,1)=angles(s,1)+180;
                    else
                        angles(s,1)=angles(s,1)-180;
                    end
%                 else 
%                     angles(s,2)=90-angles(s,2);
                end
                
                if angles(s,2)>70 
                    lab_vert=[lab_vert; s];
                    mfeFINALE(end)=[];
                end
                %                 theta=angles(s,1)*pi/180;
                %                 phi=angles(s,2)*pi/180;
                %                 out_aus=[cos(theta)*sin(phi) sin(theta)*sin(phi) cos(phi)];
                %                 quiver3(p(1),p(2),p(3),out_aus(1),out_aus(2),out_aus(3))
            end

            segmentation={};

            for k=1:numel(segment)
                if numel(find(k==lab_vert))==0
                    segmentation=[segmentation; segment{k,1}];
                end
            end
            
            Final(lab_vert,:)=[];
            angles(lab_vert,:)=[];
            angles_orig(lab_vert,:)=[];

            segment=segmentation;
            if size(segment,1)>1
                % Check oversegment
                dist=zeros(size(segment,1),size(segment,1));
                for k=1:numel(segment)
                    v1=Final(k,1:3);
                    aus1=segment{k,1};
                    p1=aus1(1,:);
                    for s=k+1:numel(segment)
                        aus2=segment{s,1};
                        v2=Final(s,1:3);
                        p2=aus2(1,:);
                        [I,D] = knnsearch(aus1,aus2);
                        %                     dist(k,s)=abs(sum(v1(1:3).*(p1-p2)))+min(D);
                        dist(k,s)=norm(cross(v1,v2))+min(D);
                        dist(s,k)=dist(k,s);
                    end
                end

                y = squareform(dist);
                Z = linkage(y,'complete');
%                             dendrogram(Z)

                % Prima era cut=mean(mean(dist))*0.1; e si può ripristinare!!!
                if mean(mean(dist))<10
                    cut=mean(mean(dist))*0.1;
                else
                    cut=mean(mean(dist))*0.01;
                end


                %             dendrogram(Z,'ColorThreshold',cut)
                T = cluster(Z,'cutoff',cut,'Criterion','distance');

                segm=cell(max(T),1);

                for s=1:max(T)
                    segm{s,1}=find(T==s);
                end

                segmentation={};
                out_fin=[];
                angle_fin=[];
                angle_fin_orig=[];
                for s=1:size(segm,1)
                    aus=segm{s,1};
                    aus_xyz=segment{aus(1),1};
                    for k=2:numel(aus)
                        aus_xyz=[aus_xyz;segment{aus(k),1}];
                    end
                    if numel(aus)>1
                        [output, mfe, thetaf, phif]=SearchPlane(aus_xyz);
                        aus2=output;
                        thetaf=thetaf*180/pi;
                        phif=phif*180/pi;
                        aus4=[thetaf phif];
                        if phif>90
                            phif=180-phif;
                            if  thetaf<180
                                thetaf=thetaf+180;
                            else
                                thetaf=thetaf-180;
                            end
%                         else
%                             phif=90-phif;
                        end
                        aus3=[thetaf phif];
                        %                     figure
                        %                     hold on
                        %                     scatter3(aus_xyz(:,1),aus_xyz(:,2),aus_xyz(:,3),'.', 'MarkerEdgeColor', [rand rand rand])
                    else
                        aus2=Final(aus(1),:);
                        aus3=angles(aus(1),:);
                        aus4=angles_orig(aus(1),:);
                    end
                    segmentation=[segmentation; aus_xyz];
                    out_fin=[out_fin; aus2];
                    angle_fin=[angle_fin; aus3];
                    angle_fin_orig=[angle_fin_orig; aus4];
                end
% 
%                                 figure
%                                 hold on
%                 
%                                 for s=1:size(segmentation,1)
%                                     aus=segmentation{s,1};
%                                     scatter3(aus(:,1),aus(:,2),aus(:,3),'.', 'MarkerEdgeColor', [rand rand rand])
%                                 end
%                                 axis equal

                segment=segmentation;
            end

            polygons={};
            lab_pol=[];
            panel_tresh=1.7;
            for s=1:size(segment,1)
                aus=segment{s,1};
                if size(aus,1)>20
                    line=zeros(100,2);
                    
                    n=out_fin(s,1:3);
                    R = 1; % depends on your data
                    density=zeros(size(aus,1),1);
                    for idx=1:size(aus,1)
                        Distances = sqrt( sum( (aus-aus(idx,:)).^2 ,2) );
                        Ninside   = length( find(Distances<=R) );
                        density(idx) = Ninside/(4*pi*R.^3/3);
                    end
                    %         s
                    density=mean(density);
                    tresh=density/2;
                    if tresh <10
                        tresh2=0.2;
                    else
                        tresh2=0.15;
                    end
                    epsilon=density*10^(-2)+0.1;
                    if abs(n(3))<0.99
                        z=[0 0 1];
                        R = rot_mat(z',n');
                        aus=aus*R;
                    end
                    Mbb = minBoundingBox(aus(:,1:2)');
                    
                    %         figure
                    %         axis equal
                    %         hold on
                    %         scatter(aus(:,1), aus(:,2),'.k');
                    %         plot(Mbb(1,[1:end 1]),Mbb(2,[1:end 1]),'r');
                    
                    Mbb=Mbb';
                    dist1=norm(Mbb(1,:)-Mbb(2,:));
                    dist2=norm(Mbb(3,:)-Mbb(2,:));
                    coefficients1 = polyfit([Mbb(1,1), Mbb(2,1)], [Mbb(1,2), Mbb(2,2)], 1);
                    coefficients2 = polyfit([Mbb(3,1), Mbb(2,1)], [Mbb(3,2), Mbb(2,2)], 1);
                    dir1=coefficients1(1);
                    dir2=coefficients2(1);
                    
                    if  dist1> panel_tresh && dist2>panel_tresh %&& density>=3
                        z=mean(aus(:,3));
                        k=boundary(aus(:,1),aus(:,2),0.8);
                        result = DouglasPeuckerB(aus(k,:),epsilon);
                        
                        %                         result=result(1:end-1,:);
                        [I,D] = knnsearch(result(:,1:2),Mbb);
                        I=I(1);
                        if I>1 % primo punto "coincide" con un estremo della mBB, riordinamento, ma mantengo antiorario
                            new_res=zeros(size(result));
                            new_res(1:size(result,1)-(I-1),:)=result(I:end,:);
                            new_res(size(result,1)-(I-2):end,:)=result(1:I-1,:);
                            result=new_res;
                        end
                        border=aus(k,1:2);
                        [ind,D] = knnsearch(border,result(I,1:2));
                        if ind>1 % primo punto "coincide" con un estremo della mBB, riordinamento, ma mantengo antiorario
                            new_border=zeros(size(border));
                            new_border(1:size(border,1)-(ind-1),:)=border(ind:end,:);
                            new_border(size(border,1)-(ind-2):end,:)=border(1:ind-1,:);
                            border=new_border;
                        end
                        
                        if norm(result(1,:)-result(end,:))>10^(-1)
                            result=[result; result(1,:)];
                        end
%                         
%                         figure(s)
%                         hold on
%                         axis equal
%                         scatter(aus(:,1),aus(:,2),'.', 'MarkerEdgeColor', [rand rand rand]); %scatter(border(ind_aus,1),border(ind_aus,2),'*r')
%                         scatter(aus(k,1),aus(k,2),'*r');
%                         plot(result(:,1), result(:,2))

%                         figure
%                         hold on
%                         axis equal
%                         scatter(aus(:,1),aus(:,2),'.k'); %scatter(border(ind_aus,1),border(ind_aus,2),'*r')
%                         scatter(aus(k,1),aus(k,2),'*b');
%                      %   plot(result(:,1), result(:,2))
%                         
%                          figure
%                         hold on
%                         axis equal
%                         scatter(aus(:,1),aus(:,2),'.', 'MarkerEdgeColor', [rand rand rand]); 
%                        % scatter(border(ind_aus,1),border(ind_aus,2),'*r')
%                        scatter(aus(k,1),aus(k,2),'*b');
%                         plot(result(:,1), result(:,2),'*r')
%                         plot(result(:,1), result(:,2))
                        
                        res_fin=[];
                        ind=0;
                        ind2=0;
                        indind=[];
                        for t=1:size(result,1)
                            if numel(find(ind==t))==0
                                if t<size(result,1)
                                    if t==1
                                        Pin=result(end-1,:);
                                        P0=result(t,:);
                                        P1=result(t+1,:);
                                    else
                                        Pin=result(t-1,:);
                                        P0=result(t,:);
                                        P1=result(t+1,:);
                                    end
                                else
                                    Pin=result(t-1,:);
                                    P0=result(t,:);
                                    P1=result(1,:);
                                end
                                coefficients = polyfit([P0(1), P1(1)], [P0(2), P1(2)], 1);
                                if abs(dir1-coefficients(1))<3*10^(-1) || abs(dir2-coefficients(1))<3*10^(-1)
                                    
                                    line1(:,1)=linspace(min(aus(:,1)),max(aus(:,1)),100);
                                    line1(:,2)=coefficients(1)*line1(:,1)+coefficients(2);
                                    
                                    [ind,D] = knnsearch(border,line1);
                                    ind=ind(D<0.3);
                                    ind=unique(ind);
                                    [I,D] = knnsearch(border,P0(1:2));
                                    ind_aus=I;
                                    ind=ind(find(ind==I):end);
                                    control=1;
                                    for p=1:numel(ind)-1
                                        if (ind(p+1)-ind(p))>20
                                            %  (ind(p+1)-ind(p))
                                            control=0;
                                        end
                                        if (ind(p+1)-ind(p))<=20 && control
                                            ind_aus=[ind_aus; ind(p+1)];
                                        end
                                    end
                                    
                                    %                         figure(t)
                                    %                         hold on
                                    %                         axis equal
                                    %                         scatter(aus(:,1),aus(:,2),'.', 'MarkerEdgeColor', [rand rand rand]); scatter(border(ind_aus,1),border(ind_aus,2),'*r')
                                    %                         scatter(result(:,1), result(:,2),'*b')
                                    
                                    if size(res_fin,1)>0
                                        [I,D] = knnsearch(res_fin,result(t,:));
                                    else
                                        D=0;
                                    end
                                    if D>10^(-1)
                                        if numel(find(indind>t))>0
                                            res_fin=[res_fin(find(indind<t),:);result(t,:); res_fin(find(indind>t),:)];
                                            indind=[indind(find(indind<t));t; indind(find(indind>t))];
                                        else
                                            res_fin=[res_fin;result(t,:)];
                                            indind=[indind;t];
                                        end
                                    end
                                    [I1,D] = knnsearch(result(:,1:2),border(ind_aus(1),:));
                                    [I2,D] = knnsearch(result(:,1:2),border(ind_aus(end),:));
                                    if numel(find(indind>I2))>0
                                        res_fin=[res_fin(find(indind<I2),:);result(I2,:); res_fin(find(indind>I2),:)];
                                        indind=[indind(find(indind<I2));I2; indind(find(indind>I2))];
                                    else
                                        res_fin=[res_fin;result(I2,:)];
                                        indind=[indind;I2];
                                    end
                                    ind=I1:I2-1;
                                    ind2=I1:I2;
                                else
                                    coefficients2 = polyfit([P0(1), Pin(1)], [P0(2), Pin(2)], 1);
                                    %                                     theta=atan((coefficients(1)-coefficients2(1))/(1+coefficients(1)*coefficients2(1)));
                                    %                                     theta=2*theta+pi;
                                    %                                     theta=2*pi-theta;
                                    %                                     theta*180/pi
                                    v1=P1-P0;
                                    v2=Pin-P0;
                                    dot = v1(1)*v2(1) + v1(2)*v2(2);
                                    det = v1(1)*v2(2) - v1(2)*v2(1);
                                    theta = atan2(det, dot);
                                    if theta<0
                                        theta = 2*pi+theta;
                                    end
                                    if (theta>=pi/6 && theta<3*pi/2) || theta==0
                                        %                                         if numel(find(ind2==t))==0
                                        if numel(find(indind>t))>0
                                            res_fin=[res_fin(find(indind<t),:); result(t,:); res_fin(find(indind>t),:)];
                                            indind=[indind(find(indind<t)); t; indind(find(indind>t))];
                                        else
                                            res_fin=[res_fin;result(t,:)];
                                            indind=[indind;t];
                                        end
                                    else
                                        ind=[ind t+1];
                                        if numel(find(indind>t))>0
                                            res_fin=[res_fin(find(indind<t),:); result(t,:); res_fin(find(indind>t),:)];
                                            indind=[indind(find(indind<t)); t; indind(find(indind>t))];
                                        else
                                            res_fin=[res_fin;result(t,:)];
                                            indind=[indind;t];
                                        end
                                    end
                                end
                            end
                        end

                        
                        if size(res_fin,1)<=round(size(result,1)*0.4) || size(res_fin,1)<=3
                            s
                            res_fin=result;
                        end
                        if norm(res_fin(1,:)-res_fin(end,:))>10^(-1)
                            res_fin=[res_fin; res_fin(1,:)];
                        end
                        
%                         Pin=res_fin(end-2,:);
%                         P0=res_fin(end-1,:);
%                         P1=res_fin(end,:);
%                         v1=P1-P0;
%                         v2=Pin-P0;
%                         dot = v1(1)*v2(1) + v1(2)*v2(2);
%                         det = v1(1)*v2(2) - v1(2)*v2(1);
%                         theta = atan2(det, dot);
%                         if theta<0
%                             theta = 2*pi+theta;
%                         end
%                         if (theta<pi/4 || theta>3*pi/2)
%                             res_fin(end-1,:)=[];
%                         end
%
% 

%                         figure
%                         hold on
%                         axis equal
%                         scatter(aus(:,1),aus(:,2),'.k'); %scatter(border(ind_aus,1),border(ind_aus,2),'*r')
%                          scatter(aus(k,1),aus(k,2),'*b');
%                         plot(res_fin(:,1), res_fin(:,2),'LineWidth',2)
%                         plot(res_fin(:,1), res_fin(:,2),'*r')


                        
                        % elimino vertici duplicati
                        aus_aus=res_fin(1,:);
                        res_fin_aus=aus_aus;
                        for t=2:size(res_fin,1)
                            if norm(res_fin(t,:)-aus_aus)>10^(-2)
                                res_fin_aus=[res_fin_aus; res_fin(t,:)];
                            end  
                            aus_aus=res_fin(t,:);
                        end
                        res_fin=res_fin_aus;
% 
%                         figure(s)
%                         hold on
%                         axis equal
%                         scatter(aus(:,1),aus(:,2),'.', 'MarkerEdgeColor', [rand rand rand]); %scatter(border(ind_aus,1),border(ind_aus,2),'*r')
%                         plot(res_fin(:,1), res_fin(:,2))
%                         plot(res_fin(:,1), res_fin(:,2),'*r')

                        res_fin(end,:)=[];

                        if abs(n(3))<0.99
                            result=res_fin*R';
                        else
                            result=res_fin;
                        end
                        aus=segment{s,1};
%                         figure(s)
%                         hold on
%                         axis equal
%                         scatter3(aus(:,1),aus(:,2),aus(:,3),'.', 'MarkerEdgeColor', [rand rand rand]); scatter3(aus(k,1),aus(k,2),aus(k,3),'*r')
%                         scatter3(result(:,1), result(:,2),result(:,3),'*b')
                        %            index=[index; s];

                        % potrebbe essere un buon controllo il numero di
                        % vertici trovati insieme alla densità?
                        %                         dens=[dens; density];
%                         if density>=20 || size(res_fin,1)<=10
                            polygons=[polygons; result];
                            lab_pol=[lab_pol; s];
%                         end
                    end
                end
            end
            
            holes={};
            lab_holes=[];
            M = shaperead('../data/edifici-3820.shp');
            X=M(INDEX_ID+1).X;
            Y=M(INDEX_ID+1).Y;
            X=X-trX;
            Y=Y-trY;
            pol={};
            control=0;
            xy=[];
            for t=1:numel(X)
                if isnan(X(t))
                    xy(end,:)=[];
                    if control>0
                        pol=[pol; xy];
                    end
                    control=control+1;
                    xy=[];
                else
                    xy=[xy; X(t) Y(t)];
                end
            end
            
            tresh=6*10^(-1);
            for s=1:size(segment,1)
                aus=segment{s,1};
                n=out_fin(s,1:3);
                if abs(n(3))<0.99
                    z=[0 0 1];
                    R = rot_mat(z',n');
                end
                for t=1:size(polygons,1)
                    if lab_pol(t)~=s  
                        p=polygons{t,1};
%                         if abs(n(3))<0.99
%                             p=p*R;
                        %                         end
                        [I,D] = knnsearch(aus(:,1:2),p(:,1:2));
                        if numel(find(D<tresh))==numel(D)
                            lab_holes=[lab_holes; s];
                            holes=[holes; polygons{t,1}];
                        end
                    end
                end
                
                for t=1:size(pol,1)
                    p=pol{t,1};
                    [I,D] = knnsearch(aus(:,1:2),p(:,1:2));
                    if numel(find(D<tresh))==numel(D)
                        if abs(n(3))<0.99
                            aus=aus*R;
                            p=[p ones(size(p,1),1)*mean(aus(:,3))];
                            p=p*R';
                        else
                            p=[p ones(size(p,1),1)*mean(aus(:,3))];
                        end
                        lab_holes=[lab_holes; s];
                        holes=[holes; p];
                    end
                end
            end
            
            dirNames=strcat('./polygons/building',  int2str(INDEX_ID));
            [status, msg, msgID] = mkdir(dirNames);
            

            for s=1:numel(lab_pol)
                nameFile=strcat(dirNames, '/polygons', int2str(lab_pol(s)),'.off');
                fileID = fopen(nameFile,'w');
                fprintf(fileID,'%6s\n','OFF');
                if numel(find(lab_pol(s)==lab_holes))==0
                    fprintf(fileID,'%5d %5d %5d\n',size(polygons{s,1},1),1,0);
                    aus=polygons{s,1};
                    aus(:,1)=aus(:,1)+trX;
                    aus(:,2)=aus(:,2)+trY;
                    aus(:,3)=aus(:,3)+trZ;
                    fprintf(fileID,'%12.8f %12.8f %12.8f\n',aus');
                    fprintf(fileID,'%5d ',size(polygons{s,1},1));
                    for p=0:size(polygons{s,1},1)-1
                        fprintf(fileID,'%5d ',p);
                    end
                    fprintf(fileID,'\n');
                    fclose(fileID);
                else
                    ind=find(lab_pol(s)==lab_holes);
                    num=0;
                    for k=1:numel(ind)
                       num=num+size(holes{ind(k),1},1); 
                    end
                    num=size(polygons{s,1},1)+num;
                    fprintf(fileID,'%5d %5d %5d\n',num,numel(find(lab_pol(s)==lab_holes))+1,0);
                    aus=polygons{s,1};
                    aus(:,1)=aus(:,1)+trX;
                    aus(:,2)=aus(:,2)+trY;
                    aus(:,3)=aus(:,3)+trZ;
                    fprintf(fileID,'%12.8f %12.8f %12.8f\n',aus');
                    for k=1:numel(ind)
                        aus=holes{ind(k),1};
                        aus(:,1)=aus(:,1)+trX;
                        aus(:,2)=aus(:,2)+trY;
                        aus(:,3)=aus(:,3)+trZ;
                        fprintf(fileID,'%12.8f %12.8f %12.8f\n',aus');
                    end
                    fprintf(fileID,'%5d ',size(polygons{s,1},1));
                    for p=0:size(polygons{s,1},1)-1
                        fprintf(fileID,'%5d ',p);
                    end
                    fprintf(fileID,'\n');
                    num=size(polygons{s,1},1);
                    for k=1:numel(ind)
                        fprintf(fileID,'%5d ',size(holes{ind(k),1},1));
                        for p=0:size(holes{ind(k),1},1)-1
                            fprintf(fileID,'%5d ',p+num);
                        end
                        fprintf(fileID,'\n');
                        num=num+size(holes{ind(k),1},1);
                    end
                    fprintf(fileID,'\n');
                    fclose(fileID);
                end
            end

            angle_aus=angle_fin;
            angle_aus(angle_fin(:,1)>=90,1)=angle_fin(angle_fin(:,1)>=90,1)-90;
            angle_aus(angle_fin(:,1)<90,1)=angle_fin(angle_fin(:,1)<90,1)+270;
            angle_aus(:,1)=360-angle_aus(:,1);
            angle_fin=angle_aus;

            aus=zeros(numel(lab_pol),3);
            aus(:,1)=lab_pol;
            for s=1:numel(lab_pol)
                aus(s,2:3)=angle_fin(lab_pol(s),:);
            end
            
            nameFile=strcat(dirNames, '/angles.txt');
            fileID = fopen(nameFile,'w');
            fprintf(fileID,'# id Nord_azimut horiz_azimut');
            writematrix(aus,nameFile,'WriteMode','append')
            fclose(fileID);
            % % fprintf(fileID,'%12f %12f %12f\n',aus');

            aus=zeros(numel(lab_pol),3);
            aus(:,1)=lab_pol;
            for s=1:numel(lab_pol)
                aus(s,2:3)=angle_fin_orig(lab_pol(s),:);
            end

            nameFile=strcat(dirNames, '/angles_HT.txt');
            fileID = fopen(nameFile,'w');
            fprintf(fileID,'# id gamma beta');
            writematrix(aus,nameFile,'WriteMode','append')
            fclose(fileID);
            
            
            for s=1:numel(lab_pol)
                nameFile=strcat(dirNames, '/polygons', int2str(lab_pol(s)),'.xyz');
                fileID = fopen(nameFile,'w');
                aus=segment{lab_pol(s),1};
                aus(:,1)=aus(:,1)+trX;
                aus(:,2)=aus(:,2)+trY;
                aus(:,3)=aus(:,3)+trZ;
                fprintf(fileID,'%12.8f %12.8f %12.8f\n',aus');
                fprintf(fileID,'\n');
                fclose(fileID);
            end
            
            
            
            %             polygons=[polygons; PP1; PP2]; da mettere qui!!
            
        end
    end
    
end

%% METRICS

disp('mean MFE')
mean(mfeFINALE)
disp('median MFE')
median(mfeFINALE)
disp('std MFE')
std(mfeFINALE)
disp('mode MFE')
mode(mfeFINALE)

