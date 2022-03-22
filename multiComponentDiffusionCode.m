%% This code runs a radially symmetric spherical multi-component diffusion model over a decreasing temperature profile 
%%with a change in cooling rate.
%%The boundary conditions used are a constant boundary condition on the rim
%%and a symmetric boundary at the core
close all
clear all

%this sets up the starting temperature and pressure
R=8.314;
Tstart=950+273;%K
P=1e4;%1Gpa in bars
T=Tstart;

%% Here we initiate the initial condition
baseAlm=0.575;
basePyr=.303;
baseGrs=.101;

midAlm=.531;
midPyr=.35;
midGrs=.1055;
midX=2.0e-2;%1.99e-2;%in m
%% set up the  rim composition constant boundary condition. 
xRim=[0.4925
    0.39
    0.107];%concentration goes in order almandine, pyrope, grossular, spessartine.
%Spessartine is excluded from the equation based on charge neutrality and
%is equal to 1-almandine-pyrope-grossular.


totalLength=.05;%m
dx=1e-4;
x=0:dx:(totalLength-dx);%m, don't have a cell at r=0 because of the problem with 1/r 
r=abs(totalLength-x);
% x=x-dx/2;


start=[ones(size(x(1):dx:midX))*midAlm ones(size((midX+dx):dx:(totalLength-dx)))*baseAlm;
    ones(size(x(1):dx:midX))*midPyr ones(size((midX+dx):dx:(totalLength-dx)))*basePyr
    ones(size(x(1):dx:midX))*midGrs ones(size((midX+dx):dx:(totalLength-dx)))*baseGrs];

old=start;
new=old;

%% this section sets up a change in cooling rate at a given temperature to 2.61 degrees.
midT=750+273;%K, this is the temperature at which the linear cooling rates switches.
midTime=0;%Myr,
startTime=22;%Myr

rate1=(Tstart-midT)/(startTime-midTime);%The first cooling rate is defined as the linear path between Tstart at the starting time and midT at the middle time (usually 750 at 1030).  
rate2=2.61;%C/myr, second cooling rate.
tmodel=0;
t650=(startTime-midTime)*(60*60*24*365.25*1e6);%archaic variable name from when 650 was the middle temperature. 

for i=1:(800000)
    %calculate diffusivity
    D=calculateDiffusionMatrix(T,P,baseAlm,basePyr,baseGrs,(1-baseAlm-basePyr-baseGrs));
    %calculate constants for finite difference
    dt=dx^2/max(D(:))/9; %time step
    dtdx2=dt/dx^2;
    dtdx=dt/dx;
    
    new(:,1)=xRim;% constant boundary condition at rim
    %centered finite difference scheme for the middle adding in
    new(1,2:end-1)=old(1,2:end-1)+dtdx2*D(1,1)*(old(1,1:end-2)-2*old(1,2:end-1)+old(1,3:end))...
        +dtdx2*D(1,2)*(old(2,1:end-2)-2*old(2,2:end-1)+old(2,3:end))+...
        dtdx2*D(1,3)*(old(3,1:end-2)-2*old(3,2:end-1)+old(3,3:end));
    new(2,2:end-1)=old(2,2:end-1)+dtdx2*D(2,1)*(old(1,1:end-2)-2*old(1,2:end-1)+old(1,3:end))...
        +dtdx2*D(2,2)*(old(2,1:end-2)-2*old(2,2:end-1)+old(2,3:end))+...
        dtdx2*D(2,3)*(old(3,1:end-2)-2*old(3,2:end-1)+old(3,3:end));
    new(3,2:end-1)=old(3,2:end-1)+dtdx2*D(3,1)*(old(1,1:end-2)-2*old(1,2:end-1)+old(1,3:end))...
        +dtdx2*D(3,2)*(old(2,1:end-2)-2*old(2,2:end-1)+old(2,3:end))+...
        dtdx2*D(3,3)*(old(3,1:end-2)-2*old(3,2:end-1)+old(3,3:end));
    %for spherical coordinates, add 2/r*du/dr to the equation,
    %to symmetric boundary condition on center
    new(1,2:end-1)=new(1,2:end-1)+dtdx*D(1,1)*(-1*old(1,1:end-2)+old(1,3:end))./r(2:end-1)...
        +dtdx*D(1,2)*(-1*old(2,1:end-2)+old(2,3:end))./r(2:end-1)+...
        dtdx*D(1,3)*(-1*old(3,1:end-2)+old(3,3:end))./r(2:end-1);
    new(2,2:end-1)=new(2,2:end-1)+dtdx*D(2,1)*(-1*old(1,1:end-2)+old(1,3:end))./r(2:end-1)...
        +dtdx*D(2,2)*(-1*old(2,1:end-2)+old(2,3:end))./r(2:end-1)+...
        dtdx*D(2,3)*(-1*old(3,1:end-2)+old(3,3:end))./r(2:end-1);
    new(3,2:end-1)=new(3,2:end-1)+dtdx*D(3,1)*(-1*old(1,1:end-2)+old(1,3:end))./r(2:end-1)...
        +dtdx*D(3,2)*(-1*old(2,1:end-2)+old(2,3:end))./r(2:end-1)+...
        dtdx*D(3,3)*(-1*old(3,1:end-2)+old(3,3:end))./r(2:end-1);
    
    %symmetric boundary condition on middle, first the regular portion    
    new(1,end)=old(1,end)+dtdx2*D(1,1)*(old(1,end-1)-old(1,end))...
        +dtdx2*D(1,2)*(2*old(2,end-1)-2*old(2,end))+...
        dtdx2*D(1,3)*(2*old(3,end-1)-2*old(3,end));
    new(2,end)=old(2,end)+dtdx2*D(2,1)*(old(1,end-1)-old(1,end))...
        +dtdx2*D(2,2)*(old(2,end-1)-old(2,end))+...
        dtdx2*D(2,3)*(old(3,end-1)-old(3,end));
    new(3,end)=old(3,end)+dtdx2*D(3,1)*(old(1,end-1)-old(1,end))...
        +dtdx2*D(3,2)*(old(2,end-1)-old(2,end))+...
        dtdx2*D(3,3)*(old(3,end-1)-old(3,end));
    
        %for spherical coordinates, add 2/r*du/dr to the equation,
    %to symmetric boundary condition on center
    new(1,end)=new(1,end)+dtdx*D(1,1)*(old(1,end-1)-old(1,end))./r(end)...
        +dtdx*D(1,2)*(old(2,end-1)-old(2,end))./r(end)+...
        dtdx*D(1,3)*(old(3,end-1)-old(3,end))./r(end);
    new(2,end)=new(2,end)+dtdx*D(2,1)*(old(1,end-1)-old(1,end))./r(end)...
        +dtdx*D(2,2)*(old(2,end-1)-old(2,end))./r(end)+...
        dtdx*D(2,3)*(old(3,end-1)-old(3,end))./r(end);
    new(3,end)=new(3,end)+dtdx*D(3,1)*(old(1,end-1)-old(1,end))./r(end)...
        +dtdx*D(3,2)*(old(2,end-1)-old(2,end))./r(end)+...
        dtdx*D(3,3)*(old(3,end-1)-old(3,end))./r(end);

    
    
    old =new;%now set old to new for next time step
    tmodel=tmodel+dt;%step forward in time
    
    %% this is the place you can change temperature cooling rate
    if T>midT %find temperature we are at for next time step.
        rate=rate1;
        T=Tstart-rate*tmodel/(60*60*24*365.25*1e6);
    else
        rate=rate2;
        T=(midT)-rate*(tmodel-t650)/(60*60*24*365.25*1e6);
    end
    %% this is optional and can be changed
    if T<773 %assume no diffusion occurs below 500C.
        break
    end
    
end
%% Plot the diffusion profile

figure(1)
close 
figure(1)
aspect=5;
width=2.5;
tt=tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

colors=magma(8);
colors=colors(4:7,:);

width1=1.5;
start(:,1)=xRim;
nexttile;
hold on
a2=plot(x*1e2,new(2,:)*100,'LineWidth', width,'Color',colors(1,:));
a3=plot(x*1e2,start(2,:)*100,'k-.','LineWidth', width-.5);

set(gca,'fontsize', 16);
pbaspect([aspect 1 1])
box on;
grid on;
axis([0 totalLength*1e2 30 40])
set(gca,'LineWidth',1.5)
set(gca,'xticklabel',{[]})
leg=legend([a2 a3],'Final Garnet Profile','Initial Condition');
set(gca,'XColor','k','YColor','k')


ylabel('Pyrope')
s = sprintf('%cC', char(176));

nexttile;
hold on

d=  plot(x*1e2,new(1,:)*100,'LineWidth', width,'Color',colors(3,:));
e=plot(x*1e2,start(1,:)*100,'k-.','LineWidth', width-.5);

box on;
grid on;

set(gca,'fontsize', 16);
ylabel('Almandine')
pbaspect([aspect 1 1])
axis([0 totalLength*1e2 45 60])
set(gca,'LineWidth',1.5)
set(gca,'xticklabel',{[]})
set(gca,'XColor','k','YColor','k')


nexttile;
hold on

plot(x*1e2,new(3,:)*100,'LineWidth', width,'Color',colors(2,:))
plot(x*1e2,start(3,:)*100,'k-.','LineWidth', width-.5)

set(gca,'fontsize', 16);
ylabel('Grossular')
pbaspect([aspect 1 1])
box on;
grid on;
axis([0 totalLength*1e2 9 11])
set(gca,'LineWidth',1.5)
set(gca,'xticklabel',{[]})
set(gca,'XColor','k','YColor','k')


nexttile;
hold on

plot(x*1e2,(1-new(1,:)-new(2,:)-new(3,:))*100,'LineWidth', width,'Color',colors(4,:));
plot(x*1e2,(1-start(1,:)-start(2,:)-start(3,:))*100,'k-.','LineWidth', width-.5);

set(gca,'fontsize', 16);
ylabel('Spessartine')
pbaspect([aspect 1 1])
axis([0 totalLength*1e2 0 3])
set(gca,'XColor','k','YColor','k')

xlabel('Distance to rim [cm]')

box on;
grid on;
set(gca,'LineWidth',1.5)
set(gcf,'Position',[0.4400    0.1665    1.1070    0.6315]*1e3)
set(leg,'Position',[0.5829    0.8246    0.2899    0.1373])

