close all; clear all;
Re = 3200;
fil = dir('*.out'); 
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
timevy = [];
timexp = [];
timexpc = [];
timeypc = [];
timeup = [];
timevp = [];
timeyp = [];
timeRnv = [];
timeRnu = [];
timeP =[];
k = 0;
for i = 1 : length(fil)
    name = fil(i).name;
    k = k + 1;
    disp(k/length(fil)*100);
    if (name(1:5) == 'vel_x')
        timevx = [timevx str2double(name (6:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vx = cat(3,Vx,A);
    elseif (name(1:5) == 'vel_y')
        timevy = [timevy str2double(name (6:end-4))];
        A = importdata(name);
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
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vxpc = cat(3,Vxpc,A);
    elseif (name(1:7) == 'vel_Ypc')
        timeypc = [timeypc str2double(name (8:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vypc = cat(3,Vypc,A);
    elseif (name(1:6) == 'vel_Xp')
        timexp = [timexp str2double(name (7:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vxp = cat(3,Vxp,A);
    elseif (name(1:6) == 'vel_Yp')
        timeyp = [timeyp str2double(name (7:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vyp = cat(3,Vxp,A);
    elseif (name(1:6) == 'vel_up')
        timeup = [timeup str2double(name (7:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vup = cat(3,Vup,A);
    elseif (name(1:6) == 'vel_vp')
        timevp = [timevp str2double(name (7:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Vvp = cat(3,Vvp,A);
    elseif (name(1:3) == 'Rnv')
        timeRnv = [timeRnv str2double(name (4:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Rnv = cat(3,Rnv,A);
    elseif (name(1:3) == 'Rnu')
        timeRnu = [timeRnu str2double(name (4:end-4))];
        A = importdata(name);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        Rnu = cat(3,Rnu,A);
    elseif (name(1:8) == 'Pressure') 
        timeP = [timeP str2double(name (9:end-4))];
        A = importdata(name);
        A(1,1)=A(1,2); A(1,end) = A(1,end-1);
        A(end,1) = A(end,2); A(end,end) = A(end,end-1);
        j = find(name(1:end-4) == '.') ;
        B = name(1:end-4);
        B(j) = 'p';
        c = [B,'=A;'];
        eval(c)
        P = cat(3,P,A);
%         subplot(1,3,3)
%         contourf(A);
    end
    
end
[fil col] = size(P(:,:,1));
fil = fil - 2;
col = col - 2;
[timevx Ivx] = sort(timevx);
filename = ['velx' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['velx' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
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
filename = ['vely' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['vely' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
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
filename = ['Pressure' num2str(fil) '_' num2str(col) '_Re' num2str(Re) '.gif'];
filenamef = ['Pressure' num2str(fil) '_' num2str(col) '_Re' num2str(Re)];
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
close all;
%%
y = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0];
sol_Rex100 = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.21090 -0.15662 -0.1015 -0.06434 -0.04775 -0.04192 -0.03717 0];
sol_Rex400 = [1 0.7583 0.68439 0.61756 0.55892 0.29093 0.16256 0.02135 -0.11477 -0.17119 -0.32726 -0.24299 -0.14612 -0.10338 -0.09266 -0.08186 0];
sol_Rex1000 = [1 0.65928 0.57492 0.51117 0.46604 0.33304 0.18719 0.05702 -0.0608 -0.10648 -0.27805 -0.38289 -0.2973 -0.2222 -0.20196 -0.18109 0];
sol_Rex3200 = [1 0.53236 0.48296 0.46547 0.46101 0.34682 0.19791 0.07156 -0.04272 -0.86636 -0.24427 -0.34323 -0.41933 -0.37827 -0.35344 -0.32407 0];
sol_Rex10000 = [1 0.47221 0.47783 0.48070 0.47804 0.34635 0.20673 0.08344 0.0311 -0.07540 -0.23186 -0.32709 -0.38 -0.41657 -0.42537 -0.42735 0];
x = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];
sol_Rey100 = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.10091 0.09233 0];
sol_Rey400 = [0 -0.12146 -0.15663 -0.19254 -0.22847 -0.23827 -0.44993 -0.38598 0.05186 0.30174 0.30203 0.28124 0.22965 0.20920 0.19713 0.18360 0];
sol_Rey1000 = [0 -0.21388 -0.27669 -0.33714 -0.39188 -0.5155 -0.42665 -0.31966 0.02526 0.32235 0.33075 0.37095 0.32627 0.30353 0.29012 0.27485 0];
sol_Rey3200 = [0 -0.39017 -0.47425 -0.52357 -0.54053 -0.44307 -0.37401 -0.31184 0.00999 0.28188 0.29030 0.37119 0.42768 0.41906 0.40917 0.39560 0];
sol_Rey10000 = [0 -0.54302 -0.52987 -0.49099 -0.45863 -0.41496 -0.36787 -0.30719 0.00831 0.27224 0.28003 0.3507 0.41487 0.43124 0.4373 0.43983 0];

fil = fil + 2;
col = col + 2;
dx = 1 / (col - 2);
dy = 1 / (fil - 2);
yvec = zeros(fil,1);
xvec = zeros(1,col);
yvec(1) = 0; yvec(2) = yvec(1) + dy / 2;
for ii = 3:fil - 1
    yvec(ii)=yvec(ii-1) + dy;
end
yvec(end) = yvec(end - 1) + dy / 2;
xvec(1) = 0; xvec(2) = xvec(1) + dx / 2;
for ii = 3:col -1
    xvec(ii) =xvec(ii - 1) + dx;
end
xvec(end) = xvec(end - 1) + dx / 2;
%%
filename = ['compx' num2str(fil-2) '_' num2str(col-2) '_Re' num2str(Re) '.gif'];
filenamef = ['compx' num2str(fil-2) '_' num2str(col-2) '_Re' num2str(Re)];
h1 = figure(1);

figure(1)
for i = 1:length(timevx)    
    plot(yvec,Vx(:,ceil(col/2),Ivx(i)),'-r', y, eval(['sol_Rex' num2str(Re)]))
    %imagesc(flipud(Vx(:,:,I(i))));
    title(['Time = ' num2str(timevx(i))])
    grid on;
    legend('Obtained Results', 'Known Results','Location','northwest')
    xlabel('y'); ylabel('u');
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

%%
filename = ['compy' num2str(fil-2) '_' num2str(col-2) '_Re' num2str(Re) '.gif'];
filenamef = ['compy' num2str(fil-2) '_' num2str(col-2) '_Re' num2str(Re)];
h1 = figure(1);

figure(1)
for i = 1:length(timevx)    
    plot(xvec,Vy(ceil(fil/2),:,Ivx(i)),'-r', x, eval(['sol_Rey' num2str(Re)]))
    %imagesc(flipud(Vx(:,:,I(i))));
    title(['Time = ' num2str(timevx(i))])
    grid on;
    legend('Obtained Results', 'Known Results','Location','northeast')
    xlabel('x'); ylabel('v');
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
