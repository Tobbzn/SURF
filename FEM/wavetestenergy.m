
%% Wave eq. test

%Length of regions a, b, c

La=20;
Lb=2;
Lc=20;
L=La+Lb+Lc;

%number of cycles for sine pulse

n_cyc=1;

%density of regions

rhoa=10;
rhob=1;
rhoc=10;

a=La/L;
b=Lb/L;
c=Lc/L;

%timesteps etc.

rhomin=min([rhoa rhob rhoc]);

Z=[5000 20000 5000];   %impedance
Zmax=max(Z);
Cmax=sqrt(Zmax/rhomin);    %max velocity 
dx=10^-2;               %preassigned dx
dt=.99*dx/Cmax;         %stability condition for dt

tmax=1.2%n_cyc;
xmax=L;
transm1=2*Z(1)/(Z(1)+Z(2));
transm2=2*Z(2)/(Z(2)+Z(3))*transm1;
refl1=(Z(1)-Z(2))/(Z(1)+Z(2));
refl2=(Z(2)-Z(3))/(Z(2)+Z(3));


A=1;                    %amplitude of sine pulse
tsteps=ceil(tmax/dt);
xsteps=ceil(xmax/dx);
u=zeros(tsteps,xsteps);

%rho=[rhoa*ones(1,xsteps*a), rhob*ones(1,xsteps*b), rhoc*ones(1,c*xsteps)];
c=sqrt([Z(1)./rhoa*ones(1,xsteps*a), Z(2)./rhob*ones(1,xsteps*b), Z(2)./rhoc*ones(1,c*xsteps)]);

ssteps=L;
mu=zeros(1,ssteps);

for m=1:ssteps

B=3*Lb/ssteps;
    
w=(pi*c(1))/(2*B*m);      %omega

s=(c*(dt/dx)).^2;       %increment scalar

T=2*pi/w;               %pulse period

u0 = @(t) A*sin(w*t);   %pulse

%initialising 
u(1,1)=u0(0);
u(2,1)=u0(dt);

%time integration
for i=3:tsteps
    u(i,1)=u0(i*dt)*(i*dt <= n_cyc*T); 
    for j=2:xsteps-1
        u(i,j)=u(i-1,j)*(2-2*s(j)) + s(j)*(u(i-1,j+1)+u(i-1,j-1)) - u(i-2,j);
    end
    
end
% 
% % movie time
% for i=1:10:tsteps
% 
%      subplot(3,1,3)
%     
%     plot(u(i,:),'r')
%     axis([ xsteps*a (xsteps*a + xsteps*b) -3*A 3*A])
%    
%     
%     subplot(3,1,1)  
%     
%     plot(u(i,:))
%     axis([ 0 xsteps -3*A 3*A])
%    
%     
%     subplot(3,1,2)
%     plot(u(i,:),'g')
%     axis([ xsteps*(a-b) (xsteps*a + xsteps*2*b) -3*A 3*A])
%     drawnow 
%     
% end

%energyplot

mu(m)=sum(mean(u.^2));
mua(m)=sum(mean(u(:,1:round(a*xsteps)).^2));
mub(m)=sum(mean(u(:,round(a*xsteps):round((a+b)*xsteps)).^2));
muc(m)=sum(mean(u(:,round((a+b)*xsteps):end).^2));
end
figure
plot(mub)

%end

