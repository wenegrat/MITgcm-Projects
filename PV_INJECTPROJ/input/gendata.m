%function gendata
%  This function generates initial condition files for the mitgcm.
% if init_vel=1, then the front problem is used
% if init_vel=0, then resting (other than small symmetry-breaking
% random noise) initial conditions are used.

%
% FLAGS FOR IC
%%%%%%%%%%%%%%%%%%%%%%%%%%
% These flags set the initial conditions to use
%IC = 'FK08'; % Single front in center of domain, following Fox-Kemper et al. 2008a
% IC = 'UNIF'; % Uniform frontal region
%IC = 'T08'; % Following Thomas and Ferrari 2008
IC = 'T08MOD'; % Modified with a tanh vertical decay of front.

% Size of domain
nx=150;
ny = 200;
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
tAlpha=-2e-4;
TtoB = -g*tAlpha;
f0=1e-4;
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
% Uncomment this stuff for stretched grid
nzs= linspace(-1,0, nz+1);
thetas = 4;
C = (1-cosh(thetas.*nzs))./(cosh(thetas)-1) ; 
hc = 300;
zf = fliplr(H*(hc*nzs + H*C)./(hc+H));
dh=-diff(zf);
zc=(zf(1:end-1)+zf(2:end))/2;  % centered z points
plot(dh, zc, '-x');
nz=length(dh);
H=sum(dh);

[XT,YT,ZT]=ndgrid(xc,yc,zc); % This is the centered, temperature grid.
[XU,YU,ZU]=ndgrid(xf,yc,zc); % This is the U grid.
[XV,YV,ZV]=ndgrid(xc,yf,zc); % This is the V grid.
[XW,YW,ZW]=ndgrid(xc,yc,zf); % This is the W grid.
[XB,YB]=ndgrid(xc,yc); % This is the Bathymetry grid.

% Stratification Params
Dml=150;             % Mixed-layer Depth
D=Dml;

y0=yc(round(length(yc)/2));     % centered location of Front

% FRONTAL PARAMATERS
Ns = (64*f0).^2;
Ns_ML =0;
Pr = nx*dxspacing*f0./sqrt(Ns)

Lf = 10000; % Half-Width of Front

% HORIZONTAL GRADIENT OPTIONS
% ==============================================================
% ==============================================================
  Ms = -(6*f0).^2; fstring='6f';
% Ms = -(3*f0).^2;
 Ms = -(1*f0).^2; fstring ='1f';
 Ms = -(2*f0).^2; fstring = '2f';
  Ms = -(4*f0).^2; fstring='4f';
% Ms = -(1.5*f0).^2;
dtdy = Ms./(9.81*2e-4);
deltatheta = dtdy.*Lf*2
MaxgradB = deltatheta*TtoB./(2*Lf) % Consistency check

T0=15; % Background temperature -not of dynamical signficance
dthetadz=-Ns/g/tAlpha;
dthetadz_ML = -Ns_ML/g/tAlpha;

fprintf(' Ms=%15.6e , Ns=%15.6e , T0= %6.3f\n',Ms,Ns,T0);
fprintf(' deltaRho=%15.6e , dthetadz=%15.6e\n',deltatheta*1035*tAlpha,dthetadz);

%MLI Scales (Fox-Kemper et al. 2008a)
Us = Ms*Dml./f0;
Ri = Ns_ML.*f0^2./Ms.^2;
% Ri = Ms.^2;
Ls = 2*pi*Us./f0*sqrt((1+Ri)/(5/2))
Lsmin = 2*pi*Us./f0*sqrt(2/5) % Minimum length scale of fastest growing MLI
ts = sqrt(54/4)*sqrt(1+Ri)/f0/3600 % in hours

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
% Impose a strong initial and restoring mixed layer to Dml
iml = find(ZT(1,1,:) > -Dml, 1, 'last');

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
     case 'T08' 
         clear uinit
        ystructure = NaN(ny,1);
        ly = ceil(ny/2);
        Ly = (ny)*dyspacing;
        ystructure(1:ly+1) = 0.5.*(1-tanh(yc(1:ly+1)./Lf) + tanh( (yc(1:ly+1) - Ly./2)./Lf));
        ysright = 0.5.*(tanh( (yc(ly+1:end)-Ly./2)./Lf) - tanh( (yc(ly+1:end)-Ly)./Lf) - 1);
        ystructure(ly+1:end) = ysright-(ysright(1)-ystructure(ly+1));
        ystructfull = repmat(permute(repmat(ystructure, [1 nx]), [2 1]), [1 1 nz]);

        ystructfull(:,:,iml+1:end) = 0;

        theta(:,:,1:iml) = deltatheta.*ystructfull(:,:,1:iml) + dthetadz_ML.*(ZT(:,:,1:iml)-ZT(1,1,iml)) + T0;
        minT = min(min(min(theta(:,:,1:iml))));
        theta(:,:,iml+1:end) = deltatheta.*ystructfull(:,:,iml+1:end) + dthetadz.*(ZT(:,:,iml+1:end)-ZT(1,1,iml+1)) + minT;
        meanGradB = nanmean(gradient(squeeze(theta(1,:,1).*TtoB), dyspacing));
        ystrucyfull = gradient(ystructfull, dyspacing);
        
        for i=1:nx
            for j=1:ny
                uinit(i,j,:) = flipud(cumtrapz(fliplr(zc), flipud(squeeze(g*tAlpha./f0*deltatheta*ystrucyfull(i,j,:)))));
            end
        end
        
        dens=theta*rho*tAlpha+rho;
        densvi = trapz(zc, squeeze(dens), 3);
        etainit = -densvi./squeeze(dens(:,:,1));
      
        etainit = etainit-etainit(1);
%         etainit = etainit.*g.*squeeze(dens(:,:,1));%This converts to pressure...But ETAN in diagnostics is actually in m???
     case 'T08MOD' 
        clear uinit
        ystructure = NaN(ny,1);
        ly = ceil(ny/2);
        Ly = (ny)*dyspacing;
        ystructure(1:ly+1) = 0.5.*(1-tanh(yc(1:ly+1)./Lf) + tanh( (yc(1:ly+1) - Ly./2)./Lf));
        ysright = 0.5.*(tanh( (yc(ly+1:end)-Ly./2)./Lf) - tanh( (yc(ly+1:end)-Ly)./Lf) - 1);
        ystructure(ly+1:end) = ysright-(ysright(1)-ystructure(ly+1));
        ystructure = ystructure - 0.5;
        ystructfull = repmat(permute(repmat(ystructure, [1 nx]), [2 1]), [1 1 nz]);
        dv = 0.15*Dml;
        vstruct = (tanh((ZT+1.5*Dml)./dv)+1)/2;
%         vstruct = exp( (ZT+Dml)/30);
%         vstruct(ZT>-Dml) = 1;
%         vstruct = 1;
        ystructfull = ystructfull.*vstruct;
        theta(:,:,1:iml) = deltatheta.*ystructfull(:,:,1:iml) + dthetadz_ML.*(ZT(:,:,1:iml)-ZT(1,1,iml)) + T0;
%         minT = min(min(min(theta(:,:,1:iml))));
        minT = T0;
        
        theta(:,:,iml+1:end) = deltatheta.*ystructfull(:,:,iml+1:end) + dthetadz.*(ZT(:,:,iml+1:end)-ZT(1,1,iml+1)) + minT;
        theta(:,:,end-1) = theta(:,:,end); %N^2 at bottom should be zero for zero heat flux.
        meanGradB = nanmean(gradient(squeeze(theta(1,:,1).*TtoB), dyspacing));
        ystrucyfull = gradient(ystructfull, dyspacing);
        
        for i=1:nx
            for j=1:ny
                uinit(i,j,:) = flipud(cumtrapz(fliplr(zc), flipud(squeeze(g*tAlpha./f0*deltatheta*ystrucyfull(i,j,:)))));
            end
        end
        
        dens=theta*rho*tAlpha+rho;
        densvi = trapz(zc, squeeze(dens), 3);
        etainit = -densvi./squeeze(dens(:,:,1));
      
        etainit = etainit-etainit(1);
end

fprintf('      Z    , Theta: j=1 ,  j=end\n');

subplot(3,2,1)
[h,c]=contourf(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(theta(1,:,:)), 20);
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


%Spice
subplot(3,2,5)
% dens=theta*rho*tAlpha+rho;
[h,c]=contour(squeeze(YT(1,:,1:end-1))/1e3,squeeze(ZT(1,:,1:end-1)),squeeze(diff(dens(1,:,:), 1, 3)));
title('Density')
xlabel('y (km)');ylabel('z (m)')
% axis([(yc(round(length(yc)/2))-1.5*Lf)/1e3 (yc(round(length(yc)/2))+1.5*Lf)/1e3 -Dml 0])
set(gca, 'clim', [-1 1].*0.03);
colorbar

subplot(3,2,4)
[h,c]=contour(squeeze(YT(1,:,:))/1e3,squeeze(ZT(1,:,:)),squeeze(dens(1,:,:)));
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
pert=1e-3*(pert-0.5);
for i=1:length(theta(1,1,:))
  theta(:,:,i)=theta(:,:,i)+pert;
end

%%%%%% NAME FILES HERE
fid=fopen(strcat('thetaInitial_', fstring, '.bin'),'w',ieee); fwrite(fid,theta,prec); fclose(fid);
fid=fopen(strcat('uInitial_',fstring,'.bin'),'w',ieee); fwrite(fid,uinit,prec); fclose(fid);
fid=fopen(strcat('etaInitial_',fstring,'.bin'),'w',ieee); fwrite(fid,etainit,prec); fclose(fid);


% Generate forcing
% Flux
% Q=zeros([1, size(XB)]);
% Q(:,:,:) = 200;

Q=zeros( size(XB));

Q(:,:) = 025;
fid=fopen('Qnet025.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);

Q(:,:) = 200;
fid=fopen('Qnet200.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);

Q(:,:) = 100;
fid=fopen('Qnet100.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);

Q(:,:) = 50;
fid=fopen('Qnet50.forcing','w',ieee); fwrite(fid,Q,prec); fclose(fid);

delR = dh;
fid=fopen('DelR.bin','w',ieee); fwrite(fid,delR,prec); fclose(fid);



%%
% Evaluate all scalings
%Ratio of surf buoyancy flux to eddy flux:
B0 = TtoB*50./(1035*3994);
disp(['B0/<w''b''>:  ', num2str(B0*f0./(0.06*MaxgradB.^2.*Dml.^2))]);
disp(['Jf/Jd:        ', num2str(0.06.*(B0*Dml).^(1/3).*MaxgradB.^2.*Dml/(f0.^2.*B0))]);
disp(['Jf/Jd^E:      ', num2str((B0*Dml).^(1/3)./(f0*Dml))]);

%%
% Make time varying Q
times = 1:121;
Qinit = 100;
amp = 225*2*pi./365; %Linearization of seasonal cycle from Kelly Dong 2013.

% plot(times, Qinit-amp*times)

Qperiodic = NaN(121,1);
ttime= 4;
Qperiodic(1:ttime) = Qinit;
Qperiodic(ttime+1:end) = Qinit - amp*(times(ttime+1:end)-times(ttime));
Qperiodic(end) = Qinit; %need this because t=0 interpolates between end and first.
% plot(times, Qperiodic);
[nx, ny] = size(XB);
QP = permute(repmat(Qperiodic, [1 nx ny]), [2 3 1]);
fid=fopen('QPeriodic.forcing','w',ieee); fwrite(fid,QP,prec); fclose(fid);