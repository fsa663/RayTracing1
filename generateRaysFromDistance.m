function [ray,yin0]=generateRaysFromDistance(csource,En,conditions,numrays,colray)

    cc3d1=En.cc3d(:,:,1);cc3d2=En.cc3d(:,:,2);cc3d3=En.cc3d(:,:,3);
    dnp=En.dnp;
    %choose elements where dnp=0
    cc3dx=cc3d1(~dnp);cc3dy=cc3d2(~dnp);cc3dz=cc3d3(~dnp);
    B=conditions.B;
    coeff=conditions.coeff;
    %E(i).cc3d(hi,hj,:)
    %outshell(:,:,1) all Xs

    xmin=min(cc3dx)-(max(cc3dx)-min(cc3dx));

    ymin=min(cc3dy);
    ymax=max(cc3dy);
    e1=0.1*abs(ymax-ymin);
    ymin=ymin-e1;ymax=ymax+e1;

    zmin=min(cc3dz);
    zmax=max(cc3dz);
    e1=0.1*abs(zmax-zmin);
    zmin=zmin-e1;zmax=zmax+e1;

    P=[xmin ymax 0;
    xmin ymin 0;
    xmin 0 zmin;
    xmin 0 zmax];

    dd=P-csource;
    thetas=atan(dd(:,2)./dd(:,1));
    rs=sqrt(dd(:,1).^2+dd(:,2).^2);
    phis=atan(dd(:,3)./rs);

    minth=min(thetas);maxth=max(thetas);
    minph=min(phis);maxph=max(phis);

    yin0=linspace(min(ymin,zmin),max(ymax,zmax),numrays);

    ang1=linspace(minth,maxth,numrays);
    ang2=linspace(minph,maxph,numrays);
    center1stElement=[mean(cc3dx) mean(cc3dy) mean(cc3dz)]; 
    vr=csource-center1stElement;
    r=sqrt(vr*vr');
    raypwr=B*((maxth-minth)*(maxph-minph))*(1/numrays^2)*exp(-r*1e-6*coeff);

    ray=[];
    for i=1:numrays
        for j=1:numrays
            %FT3=[sin(pi/2-ang2(j))*cos(ang1(i)) sin(pi/2-ang2(j))*sin(ang1(i)) cos(pi/2-ang2(j))];
            FT3(1,1)=cos(ang2(j))*cos(ang1(i));
            FT3(1,2)=cos(ang2(j))*sin(ang1(i));
            FT3(1,3)=sin(ang2(j));

            alpha=(xmin-csource(1))/FT3(1);
            P0=csource+alpha*FT3;
            ray(end+1,:)=[P0 FT3  raypwr 0     colray    1     i   j];
            %                            (1-3) (4-6)   (7)    (8)   (9)  (10)   (11)  (12)
            % 1-3: Source, 4-6: unit vector of the ray
            %7:power, 8:pathlength,9:empty,10:reflection order ,11=i, 12=j
        end
    end

end