clear

rhof=0.997e3; %medium density
rhop=1.05e3; %particle density
Cf=1.496e3; %medium speed of sound
vis=0.893e-3; %medium viscosity
a=1e-6/2; %particle radius
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
time=linspace(0,120,N);
dt=time(2)-time(1);

y=linspace(-W/2+0.5e-6,W/2-0.5e-6,Ny);
z=linspace(-h/2+0.5e-6,h/2-0.5e-6,Nz);
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

for ii=1:length(z)
    for jj=1:length(y)
        v2y(ii,jj)=3/16*v1^2/Cf*sin(2*ky*(y(jj)+W/2))*(1-3*z(ii)^2/(h/2)^2);
        Frad(ii,jj)=-2*ky*2*pi*a^3*Eac*coeff*sin(2*ky*y(jj));
        vy(ii,jj)=Frad(ii,jj)/(6*pi*vis*a)+v2y(ii,jj);
        
      if z(ii)>-h/2+a*0
        v2z(ii,jj)=3/16*v1^2/Cf*ky*h*cos(2*ky*(y(jj)+W/2))*(z(ii)^3/(h/2)^3-z(ii)/(h/2));
        vz(ii,jj)=Fnet/(6*pi*vis*a)+v2z(ii,jj);
      else
        v2z(ii,jj)=3/16*v1^2/Cf*ky*h*cos(2*ky*(y(jj)+W/2))*(z(ii)^3/(h/2)^3-z(ii)/(h/2));
        vz(ii,jj)=Fnet/(6*pi*vis*a)+v2z(ii,jj);
      end
    end
end

figure
pcolor(YG,ZG,sqrt(v2y.^2+v2z.^2));
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

vpy=zeros(1,N);
vpz=zeros(1,N);

ppy=zeros(1,N);
ppz=zeros(1,N);

ppy(1)=W/3; %initial particle y-position
ppz(1)=-5*a; %initial particle z-position

%calculating particle trajectory using 4th-order Runge Kutta method
for k=1:length(time)-1
        vpy(k)=interp2(YG,ZG,vy,ppy(k),ppz(k),'spline');
        vpz(k)=interp2(YG,ZG,vz,ppy(k),ppz(k),'spline');
        vyK1=vpy(k);
        vzK1=vpz(k);
        vyK2=interp2(YG,ZG,vy,ppy(k)+dt/2*vyK1,ppz(k)+dt/2*vzK1,'spline');
        vzK2=interp2(YG,ZG,vz,ppy(k)+dt/2*vyK1,ppz(k)+dt/2*vzK1,'spline');
        vyK3=interp2(YG,ZG,vy,ppy(k)+dt/2*vyK2,ppz(k)+dt/2*vzK2,'spline');
        vzK3=interp2(YG,ZG,vz,ppy(k)+dt/2*vyK2,ppz(k)+dt/2*vzK2,'spline');
        vyK4=interp2(YG,ZG,vy,ppy(k)+dt*vyK3,ppz(k)+dt*vzK3,'spline');
        vzK4=interp2(YG,ZG,vz,ppy(k)+dt*vyK3,ppz(k)+dt*vzK3,'spline');
        
        py=ppy(k)+dt/6*(vyK1+2*vyK2+2*vyK3+vyK4);
        if abs(py)>=W/2-a
            ppy(k+1)=sign(py)*(W/2)-sign(py)*(abs(py)-W/2+2*a);
        else
            ppy(k+1)=py;
        end
        
        pz=ppz(k)+dt/6*(vzK1+2*vzK2+2*vzK3+vzK4);
        if abs(pz)>=h/2-a
            ppz(k+1)=sign(pz)*(h/2)-sign(pz)*(abs(pz)-h/2+2*a);
        else
            ppz(k+1)=pz;
        end
end

figure
plot(ppy,ppz,'m.');
xlim([-W/2 W/2]);
ylim([-h/2-10e-6 h/2+10e-6]);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'FontSize',24);
set(gca,'LineWidth',1);
xticks([-W/2 0 W/2]);
xticklabels({'-W/2','0','W/2'});
yticks([-h/2 0 h/2]);
yticklabels({'-h/2','0','h/2'});