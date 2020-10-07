function [E]=preprocessElements(E)

[~,n]=size(E);

for i=1:n
    E(i).axis=E(i).axis./sqrt(E(i).axis*E(i).axis');

    if strcmp(E(i).type,'sphericalSurface')
        R=E(i).radius;
        A=E(i).aperture;
        bm=sqrt(R^2-(A/2)^2);
        E(i).spherecenter=E(i).center-E(i).axis*bm;

        d=E(i).axis*E(i).center';
        E(i).backplane=[E(i).axis -d];
        a=d-R;
        % Generate 2D points for plotting later
        for ti=0:5:360
            P=[E(i).spherecenter(1)+R*cos(ti*pi/180) E(i).spherecenter(2)+R*sin(ti*pi/180) 0];
            n2=E(i).backplane(1:end-1);d=E(i).backplane(end);
            Ptst=E(i).spherecenter;
            s1=(n2*Ptst'+d)>=0;s2=(n2*P'+d)>=0;
            if (s1~=s2);E(i).cc(end+1,:)=P;end;
        end 
       % Generate 3D points for plotting later
        ang1=0:5:360;ang2=0:5:180;
        for hi=1:length(ang1)
            for hj=1:length(ang2)
                ti=ang1(hi);tj=ang2(hj);
                P(1,1)=E(i).spherecenter(1)+R*sin(tj*pi/180)*cos(ti*pi/180);
                P(1,2)=E(i).spherecenter(2)+R*sin(tj*pi/180)*sin(ti*pi/180);
                P(1,3)=E(i).spherecenter(3)+R*cos(tj*pi/180);
                n2=E(i).backplane(1:end-1);d=E(i).backplane(end);
                Ptst=E(i).spherecenter;
                s1=(n2*Ptst'+d)>=0;s2=(n2*P'+d)>=0;
                if (s1~=s2);
                    E(i).cc3d(hi,hj,:)=P;E(i).dnp(hi,hj)=0;
                    %dnp is a flag array to indicate the points that should not be plotted
                else
                    E(i).cc3d(hi,hj,:)=[0 0 0];E(i).dnp(hi,hj)=1;
                end;
            end 
        end
    end

    if strcmp(E(i).type,'circularPlanarSurface')
        E(i).axis=E(i).axis./sqrt(E(i).axis*E(i).axis');
        ax=E(i).axis;
        f=[ax(2) -ax(1) 0];
        for ti=-10:10
            E(i).cc(end+1,:)=E(i).center+ti*f/10*E(i).aperture*0.5;
        end
    end 

    if (length(E(i).indexout)==0);E(i).indexout=1;end;
    if (length(E(i).indexcenter)==0);E(i).indexcenter=1;end;
    if length(E(i).isBoundary)==0;E(i).isBoundary=0;end 

end


end