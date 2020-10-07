clear all;clc; close all
% Run the file containing the configurations 
run('ex1.m'); 
numrays=10; %number of rays per axis 
resR=10; % factor of interpolation

rtplot=0; %2D raytracing plot: For fast ray tracing choose numrays<=10
spplot=0;%spot plot: For accurate spot shape simulation choose numrays>30
plot3dr=1;%3D raytracing plot, but very slow

[E]=preprocessElements(E);

conditions.B=1;% Total optical power of the target in W (Assumedn to be totally diffusive)
conditions.coeff=0.1;%Atmospheric Attenuation Coeffecient per km
R0_in_m=1000;% Range from target to first lens in meters

w=0*pi/180;%Horizontal: += ray is going right
w2=0*pi/180;%vertical: +=ray is going up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R0=R0_in_m*1000;%convert to mm
%[x y z] coordinates of the source
csource=-R0*[sin(pi/2-w2)*cos(w) sin(pi/2-w2)*sin(w) cos(pi/2-w2)];
[ray,yin0]=generateRaysFromDistance(csource ,E(7),conditions,numrays,1);
%{
w=0*pi/180;%Horizontal: += ray is going right
w2=5*pi/180;%vertical: +=ray is going up
csource=-R0*[sin(pi/2-w2)*cos(w) sin(pi/2-w2)*sin(w) cos(pi/2-w2)];
[ray2,yin0]=generateRaysFromDistance(csource ,E(7),conditions,numrays,2);
ray=[ray;ray2];
%}
clrz='rgbykmc';

% RayM is a stack for all unporcessed rays. It will be initialized to the rays coming from the source
RayM=ray; 
while length(RayM>1)
    %Take the ray at the top of the stack
    rayp=RayM(1,:);
    % Find the next intersection, and the resultant rays
    [rayReflected,rayRefracted,raycmp]=extendRay(E,rayp);
    %Collect all the completed rays
    raycmpM(end+1,:)=raycmp; 
    %Need a threshold condition on the reflected rays, otherwise the code will follow very weak rays resonating between two surfaces
    if length(rayReflected)>1; if rayReflected(7)>0.01;RayM(end+1,:)=rayReflected;end;end;
    if length(rayRefracted)>1; RayM(end+1,:)=rayRefracted;end;
    % Reduce the stack by removing the first ray (which was just processed)
    RayM=RayM(2:end,:);
end

if rtplot
    figure;
    [~,num4]=size(E);
    for i=1:num4
        if strcmp(E(i).type,'sphericalSurface')
            xx0=E(i).cc3d(:,:,1);zz0=E(i).cc3d(:,:,3);
            dnp=E(i).dnp;
            xx=xx0(~dnp);zz=zz0(~dnp);
            scatter(xx(:),zz(:));hold on;
        end
    end
    [num5,~]=size(raycmpM);
    for i=1:num5
        plot(raycmpM(i,[1 13]),raycmpM(i,[3 15]),clrz(raycmpM(i,9)),'linewidth',5*raycmpM(i,7));    
    end
    axis equal;
    hold off;title('Ray Tracing'),xlabel('x'),ylabel('z');
end                


raySens=raySensPlane=[];
%y2 raysSensor
%y2bf rays sensorplane
[a1,~]=size(raycmpM);
%Collect the completed rays falling on the sensor plane
for j=1:a1
    if (raycmpM(j,13)==E(4).center(1))
          raySensPlane(end+1,:)=raycmpM(j,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
[raySensPlaneInterp]=interpolateRays(raySensPlane,numrays,resR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

[a7,~]=size(raySensPlaneInterp);

cnd1=(raySensPlaneInterp(:,13)==E(4).center(1));
cnd2=(abs(raySensPlaneInterp(:,14))<SensLengthY/2);
cnd3=(abs(raySensPlaneInterp(:,15))<SensLengthZ/2);
cnd=logical(cnd1.*cnd2.*cnd3);

raySens=raySensPlaneInterp(cnd,:);
%{
figure;rectangle('position', [-SensLengthY/2 -SensLengthZ/2 SensLengthY SensLengthZ],'facecolor','g'); hold on;
scatter(raySensPlane(:,14),raySensPlane(:,15)),
title('Sensor Plane: Intersection of rays (Original) '), xlabel('y'),ylabel('z');hold off;
%}
%{
figure;rectangle('position', [-SensLengthY/2 -SensLengthZ/2 SensLengthY SensLengthZ],'facecolor','g');hold on;
scatter(raySensPlaneInterp(:,14),raySensPlaneInterp(:,15)),
title('Sensor Plane: Intersection of rays (Interpolated) '), xlabel('y'),ylabel('z');hold off;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy0=yin0(2)-yin0(1);
pixelSize=0.05;%image resolution
numpixels=ceil(((max(SensLengthY,SensLengthZ)/2)+dy0)/pixelSize)*2+1;
mid=round(numpixels/2-1);
xt=((1:numpixels)-mid)*pixelSize;
IM=IMdots=bx=zeros(numpixels,numpixels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a1,~]=size(raySens);
for j=1:a1
    in0x=round(mid+(raySens(j,14)/pixelSize));in0y=round(mid+(raySens(j,15)/pixelSize));
    IMdots(in0x,in0y)=IMdots(in0x,in0y)+raySens(j,7);
end

bxhw=floor(0.5*dy0/(resR*pixelSize));% box half width
bx(mid-bxhw:mid+bxhw,mid-bxhw:mid+bxhw)=ones(2*bxhw+1,2*bxhw+1);
bx=bx/sum(sum(bx));
IM=conv2(IMdots,bx,'same');


[xx,yy]=meshgrid(xt,xt');


if spplot
      figure;
      imagesc(xt,xt,IM');hold on;axis equal;
      rectangle('position', [-SensLengthY/2 -SensLengthZ/2 SensLengthY SensLengthZ]);
      xlabel('y (mm)');ylabel('z (mm)');
end

if plot3dr
    drawelements=1;Display3dLayout(raycmpM,E,drawelements);
end

