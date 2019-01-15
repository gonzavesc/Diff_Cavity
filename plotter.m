close all; clear all;
Pr = 100;
Ray = 10;
Re = 3200;
fil = dir('Results/*.out'); 
timevx = [];
k = 1;

Vx=[];
Vxp = [];
Vxpc = [];
Vypc = [];
Vyp = [];
Vup = [];
Vvp = [];
Vy=[];
Rnv = [];
Rnu = [];
P=[];
T = [];
timevy = [];
timexp = [];
timexpc = [];
timeypc = [];
timeup = [];
timevp = [];
timeyp = [];
timeRnv = [];
timeRnu = [];
timeP = [];
timeT = [];
k = 0;
for i = 1 : length(fil)
    name = fil(i).name;
    k = k + 1;
    disp(k/length(fil)*100);
    if (name(1:5) == 'vel_x')
        timevx = [timevx str2double(name (6:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vx = cat(3,Vx,A);
    elseif (name(1:5) == 'vel_y')
        timevy = [timevy str2double(name (6:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vy = cat(3,Vy,A);
%         subplot(1,3,2)
%         contourf(A);
    elseif (name(1:7) == 'vel_Xpc')
        timexpc = [timexpc str2double(name (8:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vxpc = cat(3,Vxpc,A);
    elseif (name(1:7) == 'vel_Ypc')
        timeypc = [timeypc str2double(name (8:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vypc = cat(3,Vypc,A);
    elseif (name(1:6) == 'vel_Xp')
        timexp = [timexp str2double(name (7:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vxp = cat(3,Vxp,A);
    elseif (name(1:6) == 'vel_Yp')
        timeyp = [timeyp str2double(name (7:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vyp = cat(3,Vxp,A);
    elseif (name(1:6) == 'vel_up')
        timeup = [timeup str2double(name (7:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vup = cat(3,Vup,A);
    elseif (name(1:6) == 'vel_vp')
        timevp = [timevp str2double(name (7:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vvp = cat(3,Vvp,A);
    elseif (name(1:3) == 'Rnv')
        timeRnv = [timeRnv str2double(name (4:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        j = find(name(1:end-4) == '-') ;
        B(j) = 'm';
        c = [B,'=A;'];
        eval(c)
        Rnv = cat(3,Rnv,A);
    elseif (name(1:3) == 'Rnu')
        timeRnu = [timeRnu str2double(name (4:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        j = find(name(1:end-4) == '-') ;
        B(j) = 'm';
        c = [B,'=A;'];
        eval(c)
        Rnu = cat(3,Rnu,A);
    elseif (name(1:8) == 'Pressure') 
        timeP = [timeP str2double(name (9:end-4))];
        names = "Results/" + name;
        A = importdata(names);
        A(1,1)=A(1,2); A(1,end) = A(1,end-1);
        A(end,1) = A(end,2); A(end,end) = A(end,end-1);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        P = cat(3,P,A);
    elseif (name(1:4) == 'Temp')
        timeT = [timeT str2double(name (5:end-4))];
        A = importdata(names);
        A(1,1)=A(1,2); A(1,end) = A(1,end-1);
        A(end,1) = A(end,2); A(end,end) = A(end,end-1);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        T = cat(3,T,A);
%         subplot(1,3,3)
%         contourf(A);
    end
    
end
[fil col] = size(P(:,:,1));
fil = fil - 2;
col = col - 2;
[timevx Ivx] = sort(timevx);
filename = ['Results/velx' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['Results/velx' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
h1 = figure(1);

figure(1)
for i = 1:length(timevx)    
    contourf(Vx(:,:,Ivx(i)), 'Showtext', 'on')
    %imagesc(flipud(Vx(:,:,I(i))));
    title(['Time = ' num2str(timevx(i))])
    colorbar;
    drawnow
    saveas(h1, [filenamef '_' num2str(i)], 'png');
    frame = getframe(h1);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
close all;
[timevy I] = sort(timevy);
filename = ['Results/vely' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['Results/vely' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
h2 = figure(2);

figure(2)
for i = 1:length(timevy)    
    contourf(Vy(:,:,I(i)), 'Showtext', 'on')
    title(['Time = ' num2str(timevy(i))])
    colorbar;
    drawnow
    saveas(h2, [filenamef '_' num2str(i)], 'png');
    frame = getframe(h2);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

[timeP I] = sort(timeP);
filename = ['Results/Pressure' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['Results/Pressure' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
close all;
h3 = figure(3);
figure(3)
for i = 1:length(timevy)    
    contourf(P(:,:,I(i)), 'Showtext', 'on')
    title(['Time = ' num2str(timeP(i))])
    colorbar;
    drawnow
    saveas(h3, [filenamef '_' num2str(i)], 'png');
    frame = getframe(h3);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

[timeT I] = sort(timeT);
filename = ['Results/Temp' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['Results/Temp' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
close all;
h4 = figure(4);
figure(4)
for i = 1:length(timevy)    
    contourf(T(:,:,I(i)), 'Showtext', 'on')
    title(['Time = ' num2str(timeT(i))])
    colorbar;
    drawnow
    saveas(h4, [filenamef '_' num2str(i)], 'png');
    frame = getframe(h4);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end



close all;
%%
[fil col] = size(Vxp(:,:,1));
dx = 1 / (col - 2);
dy = 1 / (fil - 2);
phi = zeros(fil,col);
err = 1e-9;
rms = err * 10;
dxs = dx * dx;
dys = dy * dy;
ap = 2 * (1 / dxs + 1 / dys);
distr = ones(col) * dx;
distr(1) = dx / 2; distr(end) = dx / 2;
distu = ones(col) * dy;
distu(1) = dy / 2; distu(end) = dy / 2;

while (rms>err)
    rms = 0;
    for (i = 2:fil-1)
        for(j = 2:col-1)
            ap = 2 * (1 / (distr(j)*distr(j-1)) + 1 / (distu(i)*distu(i-1)));
            prev = phi(i,j);
            phi(i,j) = ((phi(i,j+1) + phi (i,j-1)) / (distr(j)*distr(j-1)) + (phi(i+1,j) + phi(i-1,j)) / (distu(i)*distu(i-1)) - (Vxp(i+1,j,end)-Vxp(i-1,j,end))/(distu(i) + distu(i-1)) + (Vyp(i,j+1)-Vyp(i,j-1))/(distr(j) + distr(j-1))) / ap;
            rms = max(rms, abs(phi(i,j)-prev));
        end
    end
end

h = figure(1);
g = contour(phi, 25);
saveas(h, ['streamlines' num2str(fil) '_' num2str(col) '_Re' num2str(Re)], 'png');
saveas(h, ['streamlines' num2str(fil) '_' num2str(col) '_Re' num2str(Re)], 'epsc');
close all
