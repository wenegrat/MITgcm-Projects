nx=400;
ny = 400;
nz=100;


dxspacing=2500;
dyspacing=dxspacing;
Lx=dxspacing*nx;
Ly=dyspacing*ny;

fprintf(' nx= %i , ny= %i , nz= %i ; dx=%6.1f , dy=%6.1f\n', ...
          nx,ny,nz,dxspacing,dxspacing)

%-- Params
g=9.81;
f0=2*2*pi./(86400).*sind(30);
rho=1035;
day=24*60^2;
prec='real*8';
ieee='b';

H=300;            %Max depth

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
delR = dh;
fid=fopen('DelR.bin','w',ieee); fwrite(fid,delR,prec); fclose(fid);

%% Uncomment this stuff for stretched grid
% nzs= linspace(-1,0, nz+1);
% thetas = 4;
% C = (1-cosh(thetas.*nzs))./(cosh(thetas)-1) ; 
% hc = 300;
% zf = fliplr(H*(hc*nzs + H*C)./(hc+H));
% dh=-diff(zf);

zc=(zf(1:end-1)+zf(2:end))/2;  % centered z points
plot(dh, zc, '-x');
nz=length(dh);
H=sum(dh);

[XT,YT,ZT]=ndgrid(xc,yc,zc); % This is the centered, temperature grid.
[XU,YU,ZU]=ndgrid(xf,yc,zc); % This is the U grid.
[XV,YV,ZV]=ndgrid(xc,yf,zc); % This is the V grid.
[XW,YW,ZW]=ndgrid(xc,yc,zf); % This is the W grid.
[XB,YB]=ndgrid(xc,yc); % This is the Bathymetry grid.

%-- Bathymetry
% Bathymetry is on Xc, Yc
hh=ones(size(XB));
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


%% 
% Make SSH and U Initial Conditions
Lscale = 75*1e3; %Gaussian SSH Length Scale;
SSHa = .145; % TBD amplitude.
center = xc(floor(nx./2));
% umax = 
etainit = SSHa.*repmat(exp(-1/2.*((xc-center)./Lscale).^2), [nx 1]).*repmat(exp(-1/2.*((yc-center)./Lscale).^2).', [1, nx]);
[dEtadx, dEtady] = gradient(etainit, dx(1), dy(1));

uinit = -g./f0.*dEtady;
umax = max(max(abs(uinit)))
rossby = umax./(f0.*Lscale)
rossbySc = 1.64*rossby
vinit = g./f0.*dEtadx;


figure
subplot(3,1,1)
surf(xc./1000, yc./1000, etainit); shading interp
colorbar;

subplot(3,1,2)
quiver(xc./1000, yc./1000, uinit, vinit);
axis square
uinit = repmat(uinit.', [1 1 nz]);
vinit = repmat(vinit.', [1 1 nz]);


%%%%%% NAME FILES HERE
fid=fopen('uInitial.bin','w',ieee); fwrite(fid,uinit,prec); fclose(fid);
fid=fopen('vInitial.bin','w',ieee); fwrite(fid,vinit,prec); fclose(fid);
fid=fopen('etaInitial.bin','w',ieee); fwrite(fid,etainit,prec); fclose(fid);


% Generate forcing

% Make time varying Q
times = 1:300;
tc = 30;
taufac = 7;
tauamp = .1;

taux = tauamp.*(tanh((times-tc)./taufac)+1).*0.5;
taux = taux-taux(1);


taux(end) = taux(1); %need this because t=0 interpolates between end and first.

figure
plot(times, taux);
[nx, ny] = size(XB);
TX = permute(repmat(taux.', [1 nx ny]), [2 3 1]);
fid=fopen('taux.forcing','w',ieee); fwrite(fid,TX,prec); fclose(fid);