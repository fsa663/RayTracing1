function []=Display3dLayout(raycmpM,E,drawelements)
    figure;axis equal;
    [~,num4]=size(E);
    if drawelements
        for i=1:num4
            if strcmp(E(i).type,'sphericalSurface')
                Q=E(i).cc3d;
                dnp=E(i).dnp;
                hold on;
                [w1,w2,w3]=size(Q);
                for ti=1:w1-1
                    for tj=1:w2-1
                        xx=[Q(ti,tj,1) Q(ti,tj+1,1)];
                        yy=[Q(ti,tj,2) Q(ti,tj+1,2)];
                        zz=[Q(ti,tj,3) Q(ti,tj+1,3)];
                        if ((dnp(ti,tj)~=1)&&((dnp(ti,tj+1)~=1)))
                            plot3(xx,yy,zz);
                        end
                        xx=[Q(ti,tj,1) Q(ti+1,tj,1)];
                        yy=[Q(ti,tj,2) Q(ti+1,tj,2)];
                        zz=[Q(ti,tj,3) Q(ti+1,tj,3)];
                        if ((dnp(ti,tj)~=1)&&(dnp(ti+1,tj)~=1))
                            plot3(xx,yy,zz);
                        end
                    end%end ti
                end%end tj
            end%end spherical surface
        end %for elements
    end %end if drawelements
    %not implemented yet
    %if strcmp(E(i).type,'circularPlanarSurface'),scatter(E(i).cc(:,1),E(i).cc(:,2));hold on;end
    
    %plot rays
    clrz='rgykmc';
    [a11,a12]=size(raycmpM);
    for i=1:a11
        hold on;
        raycmp=raycmpM(i,:);
        crr=clrz(mod(i,6)+1);
        plot3(raycmp([1 13]),raycmp([2 14]),raycmp([3 15]),'r','linewidth',5*raycmp(7));
    end
    xlabel('x'); ylabel('y');zlabel('z');    
end %end function
