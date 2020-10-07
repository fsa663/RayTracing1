% This file is where you configure the setup

%Here you define the boundaries of ray tracing. You need 6 planes to define in order to enclose the simulation
%The center is just a point that lies on the plane defining the boundary
%The axis is the normal describing the plane
%isBoundary indicates that this plane terminates rays and does not reflect or refract

%Units are in mm
k0=10;sc=200;dd=120+25;
i=1; 
E(i).type='plane';E(i).center=[0 -sc 0];E(i).axis=[0 1 0];E(i).isBoundary=1;
i=2;
E(i).type='plane';E(i).center=[0 +sc 0];E(i).axis=[0 1 0];E(i).isBoundary=1;
i=3;
E(i).type='plane';E(i).center=[-sc 0 0];E(i).axis=[1 0 0];E(i).isBoundary=1;
%Element #4 is always the sensor plane
i=4; 
E(i).type='plane';E(i).center=[k0+dd+25 0 0];E(i).axis=[1 0 0];E(i).isBoundary=1;
SensLengthY=15;SensLengthZ=20;
i=5;
E(i).type='plane';E(i).center=[0 0 sc];E(i).axis=[0 0 1];E(i).isBoundary=1;
i=6;
E(i).type='plane';E(i).center=[0 0 -sc];E(i).axis=[0 0 1];E(i).isBoundary=1;

%Instead of defining lenses, I define the surfaces of the lenses. 

% For a spherical surface, the axis will be the direction where the positive part of the lens surface is pointing at. 
% I treat the spherical surface as a spherical section terminated by plane, which I will call the 'end plane'
%aperture: is the diameter of the circle at the intersection of the 'end plane' and the spherical surface
%center: is the point at the center of the circle at the intersection of the 'end plane' and the spherical surface
%radius : Radius of curvature of the spherical surface
% indexcenter: refractive index of the side of the spherical surface closest to the center of the sphere defining the surface
% indexout : refractive index of the other side of the spherical surface 
%Element #7 is always the outermost optical element facing the target. This will be used to define the rays.

%%%%%%%%%%%%%%%%%%%%%%%%FIRST LENS%%%%%%%%%%%%%%%%%%%%%%%

i=7;
E(i).type='sphericalSurface';
E(i).center=[k0 0 0]; E(i).axis=[-1 0 0];
E(i).aperture=100;E(i).radius=120;
E(i).indexout=1;E(i).indexcenter=1.5;

i=8;
E(i).type='sphericalSurface';
E(i).center=[k0 0 0]; E(i).axis=[1 0 0];
E(i).aperture=100;E(i).radius=120;
E(i).indexout=1;E(i).indexcenter=1.5;
%{
psrf1=(E(7).indexcenter-E(7).indexout)/E(7).radius;
psrf2=(E(8).indexcenter-E(8).indexout)/E(8).radius;
Plens1=psrf1+psrf2-psrf1*psrf2*(E(8).center(1)-E(7).center(1))/E(8).indexout;
flens1=1/Plens1;
%}
i=9;
E(i).type='circularAperture';
E(i).center=[k0 0 0]; E(i).axis=[-1 0 0];
E(i).aperture=100;


%%%%%%%%%%%%%%%%%%%%%%%%SECOND LENS%%%%%%%%%%%%%%%%%%%%%%%

i=10;
E(i).type='sphericalSurface';
E(i).center=[dd+k0 0 0]; E(i).axis=[-1 0 0];
E(i).aperture=20;E(i).radius=25;
E(i).indexout=1;E(i).indexcenter=1.5;
i=11;
E(i).type='sphericalSurface';
E(i).center=[dd+k0 0 0]; E(i).axis=[1 0 0];
E(i).aperture=20;E(i).radius=25;
E(i).indexout=1;E(i).indexcenter=1.5;
%{
psrf1=(E(10).indexcenter-E(10).indexout)/E(10).radius;
psrf2=(E(11).indexcenter-E(11).indexout)/E(11).radius;
Plens2=psrf1+psrf2-psrf1*psrf2*(E(11).center(1)-E(10).center(1))/E(11).indexout;
flens2=1/Plens2
%}
i=12;
E(i).type='circularAperture';
E(i).center=[dd+k0 0 0]; E(i).axis=[-1 0 0];
E(i).aperture=20;