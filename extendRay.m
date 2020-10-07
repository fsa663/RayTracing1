function [rayReflected,rayRefracted,raycmp]=extendRay(E,ray)
    cnt=1;e0=1e-12;
    [~,n]=size(E);
    for i=1:n
        F=[ray(4);ray(5);ray(6)];c0=[ray(1);ray(2);ray(3)];
         %Cycle through the optical elements, check if we have an intersection with this given ray
        if strcmp(E(i).type,'sphericalSurface')
            R=E(i).radius;
            k1=c0(1)-E(i).spherecenter(1);
            k2=c0(2)-E(i).spherecenter(2);
            k3=c0(3)-E(i).spherecenter(3);
            a=F(1)^2+F(2)^2+F(3)^2;
            b=2*k1*F(1)+2*k2*F(2)+2*k3*F(3);
            c=k1^2+k2^2+k3^2-R^2;
            del=sqrt(b^2-4*a*c);
            alpha1=(-b+del)/(2*a);
            alpha2=(-b-del)/(2*a);
            n2=E(i).backplane(1:end-1);d=E(i).backplane(end);
            Ptst=E(i).spherecenter;
            
            P=c0+alpha1*F;
            s1=(n2*Ptst'+d)>=0;s2=(n2*P+d)>=0;
            cm1=(s1~=s2)&&(alpha1>1e-12)&&(alpha1==real(alpha1));
            if (cm1) ;P1=P;end

            P=c0+alpha2*F;
            s1=(n2*Ptst'+d)>=0;s2=(n2*P+d)>=0;
            cm2=(s1~=s2)&& (alpha2>1e-12)&&(alpha2==real(alpha2));
            if (cm2);P1=P;end

            cm=(cm1||cm2);
        end%%%%%end lens

        if strcmp(E(i).type,'plane')
            nv=E(i).axis;d=nv*E(i).center';
            if (nv*F==0);cm=0;else
                alpha=(d-nv*c0)/(nv*F);
                P1=c0+alpha*F;
                h1=P1-c0;r=sqrt(h1'*h1);
                cm=(alpha>0);
            end
        end%end plane


        if strcmp(E(i).type,'circularPlanarSurface')
            nv=E(i).axis;d=nv*E(i).center';
            if (nv*F==0);
                cm=0;
            else
                alpha=(d-nv*c0)/(nv*F);
                P1=c0+alpha*F;
                h1=P1-c0;r=sqrt(h1'*h1);
                cm1=(alpha>0);
            end
            rcv=sqrt(sum((P1-E(i).center').^2));
            cm2=(rcv<=0.5*E(i).aperture);
            cm=(cm1&&cm2);
        end%%%%%end circularPlanarSurface

        if strcmp(E(i).type,'circularAperture')
            nv=E(i).axis;d=nv*E(i).center';
            if (nv*F==0);
                cm=0;
            else
                alpha=(d-nv*c0)/(nv*F);
                P1=c0+alpha*F;
                h1=P1-c0;r=sqrt(h1'*h1);
                cm1=(alpha>0);
            end
            cm=cm1;
        end%%%%%end circularaperture
        
       %If there is any possible intersection of the element with the ray,log the result
        if (cm)
            v1=P1-c0;
            %find the distance of the intersection from the source of the ray and log it
            r=sqrt(v1'*v1);
            possIntrsc(cnt,:)=[r c0'   P1' i];
            %                                  1  2-4   5-7  8 
            cnt=cnt+1;
        end
    end%%%%%%%%%%%%%%%%%%%%%%%end for

    if cnt>1
        %find the intersection closest to the source
        [~,in]=min(possIntrsc(:,1));
        crctIntrsc=possIntrsc(in,2:end);
        %find the ID of the correct element that the ray intersected
        elsrc=crctIntrsc(7);
        isBoundary=E(elsrc).isBoundary;
        c0=(crctIntrsc(1:3))';P1=(crctIntrsc(4:6))';
        F=P1-c0;
        dir=sign(E(elsrc).axis*F);
        if (dir==1); 
            nf1=E(elsrc).indexcenter;nf2=E(elsrc).indexout;
        else
            nf1=E(elsrc).indexout;nf2=E(elsrc).indexcenter; 
        end
        %dp=(F'*nv);%has to point at different diretions
        %nv=-sign(dp)*nv;

        if strcmp(E(elsrc).type,'sphericalSurface')
            nv=(P1-E(elsrc).spherecenter');nv=nv/sqrt(nv'*nv);
        end
        if ((strcmp(E(elsrc).type,'plane')) || (strcmp(E(elsrc).type,'circularPlanarSurface')))
            nv=E(elsrc).axis';nv=nv/sqrt(nv'*nv);
        end
        if (strcmp(E(elsrc).type,'circularAperture'))
            nv=E(elsrc).axis';nv=nv/sqrt(nv'*nv);
            %override isbounary for intersections outside the aperture (when killray==1)
            rcv=sqrt(sum((P1-E(elsrc).center').^2));
            isBoundary=(rcv>=0.5*E(elsrc).aperture);
        end
        %Initialize Total Internal Reflection flag to 0
        TIR=0;
        Fu0=F/sqrt(F'*F);
        th1=acos(abs(-nv'*Fu0));
        th2=asin(abs(nf1*sin(th1)/nf2));
        if abs(nf1*sin(th1)/nf2)>1, TIR=1;F22=[0;0;0];end;
        Rs=((nf1*cos(th1)-nf2*cos(th2))/(nf1*cos(th1)+nf2*cos(th2)))^2;
        Rp=((nf1*cos(th2)-nf2*cos(th1))/(nf1*cos(th2)+nf2*cos(th1)))^2;
        % For the time being, I do not factor polarization into calculation. Therefore the coeffecient of reflection is the average for TE and TM
        cR=0.5*(Rs+Rp);
        cT=1-cR;

        Fn=(F'*nv)*nv;
        F21=F-2*Fn;
        ampF=sqrt(F21'*F21); if ampF==0; ampF=1;end;
        F21=F21/ampF;
        %align nv to f
        F=F/sqrt(F'*F);
        nv2=sign(F'*nv)*nv;%aligned nv
        mv=F-(F'*nv2)*nv2;%vector
        mvamp=sqrt(mv'*mv);mvunit=mv/mvamp;
        b1=mvamp;
        thold=th1;
        th1=asin(b1);
        th2=asin(abs(nf1*sin(th1)/nf2));
        a2=cos(th2);b2=sin(th2);
        F22=a2*nv2+b2*mvunit;
        raycmp=[ray P1'];
        
        %Generate next rays (reflected and refracted)
        if (isBoundary)
            rayReflected=rayRefracted=[];
        else 
            d=sqrt(sum((P1-c0).^2));
            %reflected
            rayReflected=[P1' F21' ray(7)*cR ray(8)+d ray(9) ray(10)+1 ray(11) ray(12)];
            %rayReflected=[];
            %transmitted
            rayRefracted=[P1' F22' ray(7)*cT ray(8)+d ray(9) ray(10)      ray(11) ray(12)];
            if (TIR==1), rayRefracted=[];end;
        end
    end

end