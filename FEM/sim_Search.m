%clear all;
%close all;
format long;
addpath(genpath('../Converters'));

%Declaration of parameters:
%Polymer:
X = 1;%15*10^(-6); %Length scale.

%Meshname = 'Mediummemsphere';% I mesh-detail-rekkefølge: HQMemsphere,SMemSphere Membramewsphere
%Phys_groups = [70,71,72, 133, 134]; % < Memspheres 
%E = [1,1,1,1,1];% Memspheres
%v = [0,0,0,0,0]; % Memspheres
%rho = [1,1,1,1,1];% Memspheres

Meshname = 'QCube';
Phys_groups = 0;
E = 1; %Cube
v = 0; %Cube
rho = 1; % Cube

%Meshname = 'QMembrane'; % HQMembrane,UMembrane, %MMembrane, QMembrane
%Phys_groups =  [70,71,72]; % < Membranes
%E = 10^(0)*[1,1,1];% Membranes
%v = [0,0,0]; % Membranes
%rho = [1,1,1];% Membranes


modes = 1;%0.8:0.02:1.2; %0.2:0.05:4;%0.4:0.2:1.4;%
Energies = zeros(length(modes),length(Phys_groups)+3); %zeros(length(0.2:0.2:2),1);
InitEnergy = Energies;
simnum=0;
ispulse = 0;



output_folder = 'paraview/animation/Membrane';
steps=round(pi*1000);    % Number of time steps.
NumberOfPics = steps;
granularity = 0.02;
OLT = 0.02;


%Parameters for time integration:
T0=0;                       %start time.
%szU=size(A,1);              %dimension of our system.
tic
[p, tri, tetr, szU, K1, K2, Amod, Mmod_inv, vel,Fmat_up, Fmat_low, F_acc, lowerNodes, uz_low,upperNodes, uz_up, base_omega, dt, x_plates, y_plates] = Assembly_Membrane(Meshname,Phys_groups, E,v,rho, 1, granularity);
toc
%K1 = full(K1);
% TimeIntegrator.m-prepping:
FM_up = -Mmod_inv*sum(Fmat_up,2);
FM_acc = -Mmod_inv*F_acc;


for harmonicMode = modes%0.2:0.2:2 % Harmonic mode (2n+1)
simnum = simnum+1;
%Parameters for Paraview printing:
ExtraNameNote = 'search_atwork_';

modestring = sprintf('%d',harmonicMode);
vtktitle = [Meshname '_' ExtraNameNote modestring(1:(min(3,length(modestring)))) ];
if ispulse
    vtktitle = [vtktitle 'p' num2str(ispulse)];
end
vtktitle% = 'TESTNAME';
omega = harmonicMode*base_omega;
plateDisp = @(t) -OLT*((sin(omega*t))).*(ispulse*t<(2*pi/omega));
plateAcc = @(t) -omega^2*plateDisp(t);

F_cur = @(t) (FM_up-omega^2*FM_acc)*(plateDisp(t)+plateDisp(t-dt));

wvel = omega/(harmonicMode*pi/4);
%AnalyticSolution = @(x,t) plateDisp(t+(x-1)/wvel) - plateDisp(t-(x+3)/wvel).*(wvel*t>=(x+3)) + plateDisp(t+(x-5)/wvel).*(wvel*t>=(5-x)) - plateDisp(t-(x+7)/wvel).*(wvel*t>=(7+x));

MarkerNode = find(sum(p.^2,2)==min(sum(p.^2,2)));

TimeIntegrator;
% F_last = -Fmat_up*uz_up-Fmat_low*uz_low;
% Utemp_last = Amod\F_last;
% u0 = putDirichletBackV2(Utemp_last,lowerNodes,upperNodes,uz_low,uz_up,x_plates,y_plates);
% %U = full(u0);
% U = zeros(szU,1);
% disp('Starting time integration')
% %------------------------------------------------
% %       TIME INTEGRATION (NEWMARK)
% %------------------------------------------------
% Uz=zeros(szU/3,steps);
% n=1;
% %h = waitbar(0,'Pictures taken');
% tic
% for i=1:steps
%     
% %     if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
% %         State_to_vtk(output_folder,vtktitle,n,szU,tetr(:,1:4),p,U);
% %         n=n+1;
% %     end
%     t =T0+(i-1)*dt;
%     F_cur = -Fmat_up*(uz_up+plateDisp(t))-Fmat_low*uz_low-plateAcc(t)*F_acc;
%     Utemp_cur = K2*(Utemp_last+dt*vel-0.25*dt^2*K1*Utemp_last+0.25*dt^2*(Mmod_inv*(F_cur+F_last)));
%     utemp = putDirichletBackV2(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t),x_plates,y_plates);
%     U = utemp;
%     Uz(:,i) = U(3:3:end);
%     vel = vel+ 0.5*dt*(Mmod_inv*(F_cur+F_last))-0.5*dt*K1*(Utemp_cur+Utemp_last);
%     Utemp_last = Utemp_cur;
%     F_last = F_cur;
%     %waitbar(i/steps); 
% end
% toc
%close(h);
end
lowestmode = sprintf('%d',min(modes));
highestmode = sprintf('%d',max(modes));

EnergySaver(modes',Energies, [Meshname '_energydata_' lowestmode(1:(min(3,length(lowestmode)))) '-' highestmode(1:(min(3,length(highestmode))))])
disp('Done');