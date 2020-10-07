function [y2]=interpolateRays(y2,numrays,res)
    %Ray:
    % 1-3: initial coordinates, 4-6: angles
    %7:power, 8:pathlength,9:reserved,10:reflection order ,11=i, 12=j
    %13,14,15 terminal position
    Q=[];
    rayordmax=max(y2(:,10));

    %Rearrange all the generated rays into arrays according to their reflection order and the i,j indices 
    %of the initial ray that generated them
    for t=1:rayordmax
        for i=1:numrays
            for j=1:numrays
                cnd=(y2(:,11)==i).*(y2(:,12)==j).*(y2(:,10)==t);
                idx=find(cnd,1,'first');
                if length(idx)==0,
                        Q1(i,j,:)=zeros(1,1,15);Q1(i,j,10)=t;
                else
                        Q1(i,j,:)=y2(idx,:);
                end
            end
        end
        Q(:,:,:,t)=Q1;
    end

    u0=1:numrays;
    [ux0,uy0]=meshgrid(u0,u0');
    u1=linspace(1,numrays,res*numrays-1);
    [ux1,uy1]=meshgrid(u1,u1');
    for t=1:rayordmax
        for p=1:15
            QN(:,:,p,t)=interp2(ux0,uy0,Q(:,:,p,t),ux1,uy1);
        end
        %Lower the power by res^2, since we spread each ray into res^2 rays
        QN(:,:,7,t)=QN(:,:,7,t)/(res^2); 
    end
    
    %rearrange into the initial format
    for p=1:15
        f=QN(:,:,p,:);
        y22(:,p)=f(:);
    end
    y2=y22;
end