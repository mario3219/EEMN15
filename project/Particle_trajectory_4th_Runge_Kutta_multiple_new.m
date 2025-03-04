clc,clear

rhof=0.997e3; %medium density
rhop=1.05e3; %particle density
Cf=1.496e3; %medium speed of sound
vis=0.893e-3; %medium viscosity
a=0.5e-6/2; %particle radius
W=375e-6; %channel width
 
Kappaf=1/rhof/Cf^2; %medium compressibility
Kappap=2.49e-10; %particle compressibility
Beta=Kappap/Kappaf;
f1=1-Beta;
f2=2*(rhop-rhof)/(2*rhop+rhof);
coeff=f1/3+f2/2; %acoustic contrast factor of particle

h=133e-6; %channel height
g=9.8; 
ky=pi/W; %wave number
Eac=30; %acoustic energy density
v1=2*sqrt(Eac/rhof); %acoustic velocity of fluid particles
N=1801;
Ny=375;
Nz=133;
Np=144; %number of particles
time=linspace(0,10,N);
dt=time(2)-time(1);

y=linspace(-W/2,W/2,Ny);
z=linspace(-h/2,h/2,Nz);
v2y=zeros(Nz,Ny);
Frad=zeros(Nz,Ny);
Fydrag=zeros(Nz,Ny);
Fy=zeros(Nz,Ny);
vy=zeros(Nz,Ny);
v2z=zeros(Nz,Ny);
Fzdrag=zeros(Nz,Ny);
Fz=zeros(Nz,Ny);
vz=zeros(Nz,Ny);
[YG,ZG]=meshgrid(y,z);

Fnet=-(rhop-rhof)*g*(4/3)*pi*a^3;

nz = length(z);
ny = length(y);

for i=1:nz
    for j=1:ny
        v2y(i,j)=3/16*v1^2/Cf*sin(2*ky*(y(j)+W/2))*(1-3*z(i)^2/(h/2)^2);
        Frad(i,j)=-2*ky*2*pi*a^3*Eac*coeff*sin(2*ky*y(j));
        Fydrag(i,j)=6*pi*vis*a*v2y(i,j);
        Fy(i,j)=Frad(i,j)+Fydrag(i,j);
        vy(i,j)=Fy(i,j)/(6*pi*vis*a);
        
      if z(i)>-h/2+a
        v2z(i,j)=3/16*v1^2/Cf*ky*h*cos(2*ky*(y(j)+W/2))*(z(i)^3/(h/2)^3-z(i)/(h/2));
        Fzdrag(i,j)=6*pi*vis*a*v2z(i,j);
        Fz(i,j)=Fzdrag(i,j)+Fnet;
        vz(i,j)=Fz(i,j)/(6*pi*vis*a);
      else
        v2z(i,j)=3/16*v1^2/Cf*ky*h*cos(2*ky*(y(j)+W/2))*(z(i)^3/(h/2)^3-z(i)/(h/2));
        Fzdrag(i,j)=6*pi*vis*a*v2z(i,j);
        Fz(i,j)=Fzdrag(i,j);
        vz(i,j)=Fz(i,j)/(6*pi*vis*a);
      end
    end
end

figure
pcolor(YG,ZG,sqrt(vy.^2+vz.^2));
colormap(hot);
colorbar('vertical','FontSize',20);
shading flat
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'FontSize',24);
set(gca,'LineWidth',2);
xticks([-W/2 0 W/2]);
xticklabels({'-W/2','0','W/2'});
yticks([-h/2 0 h/2]);
yticklabels({'-h/2','0','h/2'});

vpy=zeros(Np,N);
vpz=zeros(Np,N);

ppy=zeros(Np,N);
ppz=zeros(Np,N);

%initial particle y-position
for l=1:144
    ppy(l,1)=(-W/2+5e-6)+((l-(ceil(l/12)-1)*12)-1)*((W-10e-6)/11);
end

%initial particle z-position
for n=1:144
    ppz(n,1)=h/2-2.5e-6-(ceil(n/12)-1)*(h-2*2.5e-6)/11;
end

% change interp2 to griddedInterpolant solution
[YG,ZG] = ndgrid(z,y);
vyGrid = griddedInterpolant(YG,ZG,vy,'spline','none');
vzGrid = griddedInterpolant(YG,ZG,vz,'spline','none');

%calculating particle trajectory using 4th-order Runge Kutta method
for m=1:Np
    for k=1:length(time)-1
        
        % new griddedInterpolant solution
        vpy(m,k)=vyGrid(ppz(m,k),ppy(m,k));
        vpz(m,k)=vzGrid(ppz(m,k),ppy(m,k));
        
        % old interp2 solution
        %vpy(m,k)=interp2(YG,ZG,vy,ppy(m,k),ppz(m,k),'spline');
        %vpz(m,k)=interp2(YG,ZG,vz,ppy(m,k),ppz(m,k),'spline');

        vyK1=vpy(m,k);
        vzK1=vpz(m,k);

        % new griddedInterpolant solution
        vyK2 = vyGrid(ppz(m,k)+dt/2*vzK1,ppy(m,k)+dt/2*vyK1);
        vzK2 = vzGrid(ppz(m,k)+dt/2*vzK1,ppy(m,k)+dt/2*vyK1);
        vyK3 = vyGrid(ppz(m,k)+dt/2*vzK2,ppy(m,k)+dt/2*vyK2);
        vzK3 = vzGrid(ppz(m,k)+dt/2*vzK2,ppy(m,k)+dt/2*vyK2);
        vyK4 = vyGrid(ppz(m,k)+dt*vzK3,ppy(m,k)+dt*vyK3);
        vzK4 = vzGrid(ppz(m,k)+dt*vzK3,ppy(m,k)+dt*vyK3);
        
        % old interp2 solution
        %vyK2=interp2(YG,ZG,vy,ppy(m,k)+dt/2*vyK1,ppz(m,k)+dt/2*vzK1,'spline');
        %vzK2=interp2(YG,ZG,vz,ppy(m,k)+dt/2*vyK1,ppz(m,k)+dt/2*vzK1,'spline');
        %vyK3=interp2(YG,ZG,vy,ppy(m,k)+dt/2*vyK2,ppz(m,k)+dt/2*vzK2,'spline');
        %vzK3=interp2(YG,ZG,vz,ppy(m,k)+dt/2*vyK2,ppz(m,k)+dt/2*vzK2,'spline');
        %vyK4=interp2(YG,ZG,vy,ppy(m,k)+dt*vyK3,ppz(m,k)+dt*vzK3,'spline');
        %vzK4=interp2(YG,ZG,vz,ppy(m,k)+dt*vyK3,ppz(m,k)+dt*vzK3,'spline');
        
        py=ppy(m,k)+dt/6*(vyK1+2*vyK2+2*vyK3+vyK4);
        if abs(py)>=W/2-a
            ppy(m,k+1)=sign(py)*(W/2)-sign(py)*(abs(py)-W/2+2*a);
        else
            ppy(m,k+1)=py;
        end
        
        pz=ppz(m,k)+dt/6*(vzK1+2*vzK2+2*vzK3+vzK4);
        if abs(pz)>=h/2-a
            ppz(m,k+1)=sign(pz)*(h/2)-sign(pz)*(abs(pz)-h/2+2*a);
        else
            ppz(m,k+1)=pz;
        end
    end
end


figure
plot(ppy(:,:),ppz(:,:),'b.');
xlim([-W/2 W/2]);
ylim([-h/2 h/2]);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'FontSize',24);
set(gca,'LineWidth',1);
xticks([-W/2 0 W/2]);
xticklabels({'-W/2','0','W/2'});
yticks([-h/2 0 h/2]);
yticklabels({'-h/2','0','h/2'});
