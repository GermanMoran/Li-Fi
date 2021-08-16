
function [ptx1,pty1,ptx2,pty2,ptx3,pty3,ptx4,pty4]=CuadriculaTx(tx1,ty1,tx2,ty2,tx3,ty3,tx4,ty4)
    %Variables Generales
    r=0.5;            
    rx=r;ry=r;
    l=linspace(0,2*pi);

    %Transmisor 1
    xv1=rx*cos(l)'+tx1;
    yv1=ry*sin(l)'+ty1;
    ptx1=[xv1;xv1(1)];pty1=[yv1;yv1(1)];

    %Transmisor 2
    xv2=rx*cos(l)'+tx2;
    yv2=ry*sin(l)'+ty2;
    ptx2=[xv2;xv2(1)];pty2=[yv2;yv2(1)];

    %Transmisor 3
    xv3=rx*cos(l)'+tx3;
    yv3=ry*sin(l)'+ty3;
    ptx3=[xv3;xv3(1)];pty3=[yv3;yv3(1)];

    %Transmisor 4
    xv4=rx*cos(l)'+tx4;
    yv4=ry*sin(l)'+ty4;
    ptx4=[xv4;xv4(1)];pty4=[yv4;yv4(1)];

end


