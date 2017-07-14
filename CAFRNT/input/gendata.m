%function gendata
%  This function generates initial condition files for the mitgcm.
% if init_vel=1, then the front problem is used
% if init_vel=0, then resting (other than small symmetry-breaking
% random noise) initial conditions are used.

%
% FLAGS FOR IC
%%%%%%%%%%%%%%%%%%%%%%%%%%

load('model_init_forc.mat')

IC = 'T08MOD'; % Modified with a tanh vertical decay of front.

% Size of domain
nx=150;
ny = 400;
nz=50;
init_vel=1;
dxspacing=500;
dyspacing=dxspacing;
Lx=dxspacing*nx;
Ly=dyspacing*ny;

fprintf(' nx= %i , ny= %i , nz= %i ; dx=%6.1f , dy=%6.1f\n', ...
          nx,ny,nz,dxspacing,dxspacing)

%-- Params
g=9.81;
tAlpha=-init.alpha; % XXXX
sBeta= init.beta; %XXXXXX;

f0=8.55e-5;
rho0=1024;
day=24*60^2;
prec='real*8';
ieee='b';

H=150;            %Max depth

fprintf(' Lx=%6.1f , Ly=%6.1f , H=%6.1f\n',Lx,Ly,H);

%-- Grid: x
dx=ones(1,nx);                                  % uniform resolution
dx=dx*Lx/sum(dx); 
xf=cumsum([0 dx]); % Face x points
xc=(xf(1:end-1)+xf(2:end))/2; % Centered x points

%-- Grid: y
dy=ones(1,ny);                                  % uniform resolution
dy=dy*Ly/sum(dy); 
yf=cumsum([0 dy]);  % Face y-points
yc=(yf(1:end-1)+yf(2:end))/2;  %Centered y-points
L=yc(end)-yc(1);	% this takes into account the wall of topography!!!

%-- Grid: z
dh=H/nz*ones(1,nz);
zf=-round(cumsum([0 dh])/1)*1;   % Face z points
% Uncomment this stuff for stretched grid
% nzs= linspace(-1,0, nz+1);
% thetas = 4;
% C = (1-cosh(thetas.*nzs))./(cosh(thetas)-1) ; 
% hc = 300;
% zf = fliplr(H*(hc*nzs + H*C)./(hc+H));

dh=-diff(zf);
zc=(zf(1:end-1)+zf(2:end))/2;  % centered z points
plot(dh, zc, '-x');
nz=length(dh);
H=sum(dh);

%%
[XT,YT,ZT]=ndgrid(xc,yc,zc); % This is the centered, temperature grid.
[XU,YU,ZU]=ndgrid(xf,yc,zc); % This is the U grid.
[XV,YV,ZV]=ndgrid(xc,yf,zc); % This is the V grid.
[XW,YW,ZW]=ndgrid(xc,yc,zf); % This is the W grid.
[XB,YB]=ndgrid(xc,yc); % This is the Bathymetry grid.

y0=yc(round(length(yc)/2));     % centered location of Front


% HORIZONTAL GRADIENT OPTIONS
% ==============================================================
% ==============================================================

%-- Bathymetry
% Bathymetry is on Xc, Yc
hh=ones(size(XB));
% if ~ strcmp(IC, 'T08')
% hh(:,end)=0*hh(:,end);
% end
hh=-H*hh;

figure(1); clf
subplot(221)
pcolor(XB/1e3,YB/1e3,hh)
axis equal

subplot(223)
plot(xc/1e3,dx/1e3,'.'); xlabel('X (km)');ylabel('\Delta x (km)')
title('\Delta x')

subplot(222)
plot(dy/1e3,yc/1e3,'.'); ylabel('Y (km)');xlabel('\Delta y (km)')
title('\Delta y')

subplot(224)
plot(dh,zc,'.'); ylabel('Z (m)');xlabel('\Delta z (m)')
title('\Delta z')

fid=fopen('topo_sl.bin','w',ieee); fwrite(fid,hh,prec); fclose(fid);

% Initial Temp Profile
figure(2); clf

% theta=T0+dthetadz*(ZT-ZT(1,1,end))+deltatheta*tanh((YT-y0)/Lf)/2;
theta = NaN(nx, ny, nz);
sal = theta;

%%
% Define Frontal Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear uinit
Lf = 3000;
% Lf = 2000;

ystructure = NaN(ny,1);
ly = ceil(ny/2);
Ly = (ny)*dyspacing;
ystructure(1:ly+1) = 0.5.*(1-tanh(yc(1:ly+1)./Lf) + tanh( (yc(1:ly+1) - Ly./2)./Lf));
ysright = 0.5.*(tanh( (yc(ly+1:end)-Ly./2)./Lf) - tanh( (yc(ly+1:end)-Ly)./Lf) - 1);
ystructure(ly+1:end) =( ysright-(ysright(1)-ystructure(ly+1)));
ystructure = ystructure - 0.5;
ystructfull = repmat(permute(repmat(ystructure, [1 nx]), [2 1]), [1 1 nz]);
% ystructfull = -ystructfull-0.5;
ystructfull = -ystructfull;

vstruct = (tanh((ZT+40)./7.5)+1)/2;

ystructfull = ystructfull.*vstruct;

vstructT = 10.09.*exp(0.00127.*ZT) + 12.26.*exp(0.02673.*ZT);
vstructT = vstructT./repmat(10.09.*exp(0.00127.*-6) + 12.26.*exp(0.02673.*-6), [nx ny nz]);
dTdz = 10.09.*0.00127.*exp(0.00127.*ZT) + 12.26.*0.02673.*exp(0.02673.*ZT);
vstructS = -2.309.*exp(0.009806.*ZT) + 34.7.*exp(4.679e-5.*ZT);
vstructS = vstructS./(-2.309.*exp(0.009806.*-6) + 34.7.*exp(4.679e-5.*-6));


deltatheta = 0.8*2;
deltasal = 0.25*2;
T0 = 15.8;
S0 = 33.1;

theta(:,:,:) = deltatheta.*ystructfull(:,:,:) + vstructT + T0;
sal(:,:,:) = deltasal.*-ystructfull + vstructS + S0;
% contourf(squeeze(yc(1,:,1)), squeeze(ZT(1,1,:)), squeeze(gradient(theta(1,:,:),3)).')
rho = (tAlpha*(theta-T0) + sBeta*(sal-S0)).*rho0+rho0;
b = -g*rho./rho0;
contourf(squeeze(yc(1,:,1)), squeeze(ZT(1,1,:)), squeeze(rho(1,:,:)).')

%%
[dBdy, ~, dBdz] = gradient(b, dyspacing, dxspacing, dh);

uinit = NaN(nx, ny,nz);
for i=1:nx
    for j=1:ny
        uinit(i,j,:) = flipud(cumtrapz(fliplr(zc), flipud(squeeze(-1./f0.*dBdy(i,j,:)))));
    end
end

densvi = trapz(zc, squeeze(rho), 3);
etainit = -densvi./squeeze(rho(:,:,1));


etainit = 1./squeeze(b(:,:,1)).*cumtrapz(yc, squeeze(f0.*uinit(:,:,1)),2);

etainit = etainit-etainit(1);

        
fprintf('      Z    , Theta: j=1 ,  j=end\n');

subplot(3,2,1)
% [h,c]=contourf(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(theta(1,:,:)), 20);
% axis([(yc(round(length(yc)/2))-1.5*Lf)/1e3 (yc(round(length(yc)/2))+1.5*Lf)/1e3 -Dml 0])
% colorbar
% title('Potl Temp')
xlabel('y (km)');ylabel('z (m)')

subplot(3,2,2)
[h,c]=contourf(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(theta(1,:,:)));
colorbar
title('Potl Temp')
xlabel('y (km)');ylabel('z (m)')


subplot(3,2,3)
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(rho(1,:,:)));
title('Density')
xlabel('y (km)');ylabel('z (m)')
% axis([(yc(round(length(yc)/2))-1.5*Lf)/1e3 (yc(round(length(yc)/2))+1.5*Lf)/1e3 -Dml 0])


%Spice
subplot(3,2,5)
% dens=theta*rho*tAlpha+rho;
[h,c]=contour(squeeze(YT(1,:,1:end-1))/1e3,squeeze(ZT(1,:,1:end-1)),squeeze(diff(rho(1,:,:), 1, 3)));
title('Density')
xlabel('y (km)');ylabel('z (m)')
% axis([(yc(round(length(yc)/2))-1.5*Lf)/1e3 (yc(round(length(yc)/2))+1.5*Lf)/1e3 -Dml 0])
set(gca, 'clim', [-1 1].*0.03);
colorbar

subplot(3,2,4)
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(rho(1,:,:)));
title('Density')
xlabel('y (km)');ylabel('z (m)')

subplot(3,2,6)
plot(YB(1,:)/1e3,hh(1,:))
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(uinit(1,:,:)));
title('U Initial')
xlabel('y (km)');ylabel('value')
colorbar

%-- Perturb Initial temperature
pert=rand(size(theta(:,:,1)));
pert=1e-2*(pert-0.5);
for i=1:length(theta(1,1,:))
  theta(:,:,i)=theta(:,:,i)+pert;
end

%-- Perturb Initial temperature
pert=rand(size(sal(:,:,1)));
pert=1e-2*(pert-0.5);
for i=1:length(sal(1,1,:))
  sal(:,:,i)=sal(:,:,i)+pert;
end

%%%%%% NAME FILES HERE
fid=fopen(strcat('thetaInitial','.bin'),'w',ieee); fwrite(fid,theta,prec); fclose(fid);
fid=fopen(strcat('salInitial',  '.bin'),'w',ieee); fwrite(fid,sal,prec); fclose(fid);
fid=fopen(strcat('uInitial','.bin'),'w',ieee); fwrite(fid,uinit,prec); fclose(fid);
fid=fopen(strcat('etaInitial','.bin'),'w',ieee); fwrite(fid,etainit,prec); fclose(fid);


% 
% Q=zeros( size(XB));
% 
% Q(:,:) = 025;
% fid=fopen('Qnet025.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);
% 
% Q(:,:) = 200;
% fid=fopen('Qnet200.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);
% 
% Q(:,:) = 100;
% fid=fopen('Qnet100.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);
% 
% Q(:,:) = 50;
% fid=fopen('Qnet50.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);

delR = dh;
fid=fopen('DelR.bin','w',ieee); fwrite(fid,delR,prec); fclose(fid);

%%
% 
% %%
% % Evaluate all scalings
% %Ratio of surf buoyancy flux to eddy flux:
% B0 = TtoB*50./(1035*3994);
% disp(['B0/<w''b''>:  ', num2str(B0*f0./(0.06*MaxgradB.^2.*Dml.^2))]);
% disp(['Jf/Jd:        ', num2str(0.06.*(B0*Dml).^(1/3).*MaxgradB.^2.*Dml/(f0.^2.*B0))]);
% disp(['Jf/Jd^E:      ', num2str((B0*Dml).^(1/3)./(f0*Dml))]);
% 
% %% CURRENT LW ONLY
% Make time varying Q
times = 1:121;
Qinit = 100;
amp = 225*2*pi./365; %Linearization of seasonal cycle from Kelly Dong 2013.

% plot(times, Qinit-amp*times)

Qperiodic = NaN(121,1);
Q=NaN(73,1);
tx = Q;
ty = Q;

Q(1:72) = -forc.qnet(end-71:end);
Q(end) = 0;

rampfunc = 0.5.*(1+tanh(((1:73)-12)/3));
rampfunc = rampfunc - rampfunc(1);
Q = Q.*rampfunc.';


tx(1:72) = -forc.tauy(end-71:end);
ty(1:72) = forc.taux(end-71:end);
tx(end) = 0; ty(end)=0;

tx = tx.*rampfunc.';
ty = ty.*rampfunc.';

tx = permute(repmat(tx, [1 nx ny]), [2 3 1]);
ty = permute(repmat(ty, [1 nx ny]), [2 3 1]);
Q = permute(repmat(Q, [1 nx ny]), [2 3 1]);

fid=fopen('QNet.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);
fid=fopen('tx.forcing','w',ieee); fwrite(fid,tx,prec); fclose(fid);
fid=fopen('ty.forcing','w',ieee); fwrite(fid,ty,prec); fclose(fid);

%%
ttime= 5;
Qperiodic(1:ttime) = Qinit;
Qperiodic(ttime+1:end) = Qinit - amp*(times(ttime+1:end)-times(ttime));
Qperiodic(end) = Qinit; %need this because t=0 interpolates between end and first.
% plot(times, Qperiodic);
[nx, ny] = size(XB);
QP = permute(repmat(Qperiodic, [1 nx ny]), [2 3 1]);
fid=fopen('QPeriodic.forcing','w',ieee); fwrite(fid,QP,prec); fclose(fid);
% 
% %% NEW INCLUDING SW
% 
% times = 1:121;
% Qinit =100;
% amp = 300/100; %Linearization of seasonal cycle from Kelly Dong 2013.
% 
% % plot(times, Qinit-amp*times)
% 
% Qperiodic = NaN(121,1);
% ttime= 1;
% Qperiodic(1:ttime) = Qinit;
% Qperiodic(ttime+1:end) = Qinit - amp*(times(ttime+1:end)-times(ttime));
% Qperiodic(end) = Qinit; %need this because t=0 interpolates between end and first.
% % plot(times, Qperiodic);
% 
% QP = permute(repmat(Qperiodic, [1 nx ny]), [2 3 1]);
% fid=fopen('QPeriodic.forcing','w',ieee); fwrite(fid,QP,prec); fclose(fid);
% 
% QC = -50.*ones(size(QP));
% fid=fopen('Qconstant.forcing','w',ieee); fwrite(fid,QC,prec); fclose(fid);
% 
% Qsw = -120;
% amp = -200./100
% QSWP = (Qsw + amp*(times(1:end))).';
% QSWP(end) = QSWP(1);
% QSWP = permute(repmat(QSWP, [1 nx ny]), [2 3 1]);
% 
% fid=fopen('QSW.forcing','w',ieee); fwrite(fid,QSWP,prec); fclose(fid);
% 
% plot(squeeze(QP(1,1,:)-QSWP(1,1,:)));
% hold on
% plot(squeeze(QSWP(1,1,:)));
% plot(squeeze(QP(1,1,:)), 'LineWidth', 2);
% hold off
% 
% 
% 
% %%
% cl = [-1.5 1.5];
% figure
% contour(yc./1000, zc, squeeze(theta(1,:,:)).',10,'k');
% hold on
% contour(yc./1000, zc, squeeze(uinit(1,:,:)./(abs(Us))).',12, 'LineWidth', 2);
% hold off
% set(gca, 'clim', cl);
% cb = colorbar;
% set(get(cb, 'ylabel'), 'String', '$\frac{u_o}{U}$','Rotation', 0, 'Interpreter', 'Latex', 'FontSize', 24)
% set(cb, 'Ticks', cl(1):0.5:cl(end))
% hold on
% plot(yc(1)./1000.*ones(size(zc)), zc, '>', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
% hold off
% xlabel('y (km)');
% ylabel('z (m)');
% set(gca, 'FontSize', 16);
% set(gcf, 'Color', 'w', 'Position', [ 675   473   748   501]);
