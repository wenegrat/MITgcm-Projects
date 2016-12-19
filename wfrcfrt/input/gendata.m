%function geninit()
%  This function generates initial condition files for the mitgcm.
% if init_vel=1, then the front problem is used
% if init_vel=0, then resting (other than small symmetry-breaking
% random noise) initial conditions are used.

%
% FLAGS FOR IC
IC = 'FK08'; % Single front in center of domain, following Fox-Kemper et al. 2008a
% IC = 'UNIF'; % Uniform frontal region

nx=48;
ny=nx/2+1;
ny = nx;
% ny = 24;
nz=100;
init_vel=1;
dxspacing=1000;
dyspacing=dxspacing;

Lx=dxspacing*nx;
Ly=dyspacing*ny;

fprintf(' nx= %i , ny= %i , nz= %i ; dx=%6.1f , dy=%6.1f\n', ...
          nx,ny,nz,dxspacing,dxspacing)

%-- Params
g=9.81;
tAlpha=-2e-4;
f0=1e-4;
rho=1035;
day=24*60^2;
prec='real*8';
ieee='b';

H=500;            %Max depth

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
dh=-diff(zf);
zc=(zf(1:end-1)+zf(2:end))/2;  % centered z points
nz=length(dh);
H=sum(dh);

[XT,YT,ZT]=ndgrid(xc,yc,zc); % This is the centered, temperature grid.
[XU,YU,ZU]=ndgrid(xf,yc,zc); % This is the U grid.
[XV,YV,ZV]=ndgrid(xc,yf,zc); % This is the V grid.
[XW,YW,ZW]=ndgrid(xc,yc,zf); % This is the W grid.
[XB,YB]=ndgrid(xc,yc); % This is the Bathymetry grid.

% Stratification Params
Dml=100;             % Mixed-layer Depth
D=Dml;

y0=yc(round(length(yc)/2));     % centered location of Front

Ms=(2e-4)^2 ;     %db/dy  2e-4 is about 1 K/50 km
Ms = (4*f0).^2;
% Ns=9*Ms^2/f0^2 ;  %db/dz
Ns = 16*Ms;
Ns_ML = Ms;
%MLI Scales (Fox-Kemper et al. 2008a)
Us = Ms*Dml./f0;
Ri = Ns_ML.*f0^2./Ms.^2;
Ls = 2*pi*Us./f0*sqrt((1+Ri)/(5/2))
Lf=2*1e3;          % Half-Width of Front
Lf = 1*Ls;
deltatheta=Lf*Ms/g/tAlpha;
T0=17;
dthetadz=-Ns/g/tAlpha;
dthetadz_ML = -Ns_ML/g/tAlpha;

fprintf(' Ms=%15.6e , Ns=%15.6e , T0= %6.3f\n',Ms,Ns,T0);
fprintf(' deltatheta=%15.6e , dthetadz=%15.6e\n',deltatheta,dthetadz);

%-- Bathymetry
% Bathymetry is on Xc, Yc
hh=ones(size(XB));
hh(:,end)=0*hh(:,end);
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
% Impose a strong initial and restoring mixed layer to Dml
i=1;
iml=1;

iml = find(ZT(1,1,:) > -Dml, 1, 'last');
% while ((ZT(1,1,iml)>-Dml)&(iml<length(ZT(1,1,:))))
%   iml=iml+1;
% end
switch IC
    case 'FK08'
        theta(:,:,1:iml) = T0+dthetadz_ML*(ZT(:,:,1:iml)-ZT(1,1,iml))+deltatheta*tanh(2*(YT(:,:,1:iml)-y0)/Lf)/2;
        theta(:,:,iml+1:end) = T0+dthetadz*(ZT(:,:,iml+1:end)-ZT(1,1,iml))+deltatheta*tanh(2*(YT(:,:,iml+1:end)-y0)/Lf)/2;
        % Velocity
        dthetady = deltatheta*sech(2*(YT-y0)/Lf).^2/(Lf);
        uinit = g*tAlpha/f0*dthetady.*permute(repmat(zc-zc(end), [nx 1 ny]), [1 3 2]);
    case 'UNIF'
        Mst = Ms/g/tAlpha;
        theta(:,:,1:iml) = T0+dthetadz_ML*(ZT(:,:,1:iml)-ZT(1,1,iml))+Mst.*YT(:,:,1:iml);
        theta(:,:,iml+1:end) = T0+dthetadz*(ZT(:,:,iml+1:end)-ZT(1,1,iml))+Mst.*YT(:,:,iml+1:end);
        % Velocity
        dthetady = Mst;
        uinit = g*tAlpha/f0*dthetady.*permute(repmat(zc-zc(end), [nx 1 ny]), [1 3 2]);
end

fprintf('      Z    , Theta: j=1 ,  j=end\n');
fprintf(' %10.4f %10.4f %10.4f\n', ...
[ZT(1,1,1) theta(1,1,1) theta(1,end,1) ; ...
 ZT(1,1,i) theta(1,1,i) theta(1,end,i);  ...
 ZT(1,1,end) theta(1,1,end) theta(1,end,end)]' ...
);

subplot(3,2,1)
[h,c]=contourf(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(theta(1,:,:)));
axis([(yc(round(length(yc)/2))-1.5*Lf)/1e3 (yc(round(length(yc)/2))+1.5*Lf)/1e3 -Dml 0])
colorbar
title('Potl Temp')
xlabel('y (km)');ylabel('z (m)')

subplot(3,2,2)
[h,c]=contourf(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(theta(1,:,:)));
colorbar
title('Potl Temp')
xlabel('y (km)');ylabel('z (m)')


subplot(3,2,3)
dens=theta*rho*tAlpha+rho;
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(dens(1,:,:)));
title('Density')
xlabel('y (km)');ylabel('z (m)')
axis([(yc(round(length(yc)/2))-1.5*Lf)/1e3 (yc(round(length(yc)/2))+1.5*Lf)/1e3 -Dml 0])


fid=fopen('uInitial.bin','w',ieee); fwrite(fid,uinit,prec); fclose(fid);

%Spice

x1=xc(round(length(xc)/4));     % centered location of Front #1
x2=xc(round(3*length(xc)/4));     % centered location of Front #2
spice=T0+dthetadz*(ZT-ZT(1,1,end))+deltatheta*(tanh((XT-x1)/Lf)-tanh((XT-x2)/Lf)-1)/2;

i=1;
while (i<iml)
  spice(:,:,i)=spice(:,:,iml);
  i=i+1;
end

subplot(3,2,4)
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(dens(1,:,:)));
title('Density')
xlabel('y (km)');ylabel('z (m)')

subplot(3,2,5)
[h,c]=contourf(squeeze(XT(:,1,:))/1e3,squeeze(ZT(:,1,:)),squeeze(spice(:,1,:)));
title('Spice')
xlabel('x (km)');ylabel('z (m)')

subplot(3,2,6)
plot(YB(1,:)/1e3,hh(1,:))
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(uinit(1,:,:)));
title('topography')
xlabel('x (km)');ylabel('value')

%-- Perturb Initial temperature
pert=rand(size(theta(:,:,1)));
pert=1e-5*(pert-0.5);
for i=1:length(theta(1,1,:))
  theta(:,:,i)=theta(:,:,i)+pert;
end

fid=fopen('thetaInitial.bin','w',ieee);
fwrite(fid,theta,prec); fclose(fid);

fid=fopen('spiceInitial.bin','w',ieee);
fwrite(fid,spice,prec); fclose(fid);

% Generate forcing
% Flux
% Q=zeros([1, size(XB)]);
% Q(:,:,:) = 200;

Q=zeros( size(XB));
Q(:,:) = 0;200;

fid=fopen('Qnet.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);

UW=zeros( size(XB));
UW(:,:) = .1;

fid=fopen('ZonalTau.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);
