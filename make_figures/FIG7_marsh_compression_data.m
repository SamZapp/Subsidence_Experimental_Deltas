%Reproduce Fig7

%load in data to workspace from figshare: 
load pordat.mat
load subs19_10hr_mw.mat

elev_mid=pordat(:,4)-((pordat(:,5)+pordat(:,6))./2);
n=pordat(:,16);%porosity
e=n./(1-n);%void ratio

%add in initial conditions: n=0.9,e=9,z=0 
Z=[elev_mid; 0.025];
n2=[n; 0.9];
e2=[e; 9];
%%%%%%%%%%%%%%
figure
plot(e2,Z,'o')
axis ij
xlabel('void ratio (e)')
ylabel('depth (cm)')
title('marsh compression')

%fit curve (probably log) to data
showfit exp
makevarfit
a
b
ap=a;
bp=b;
hold on

%model subsidence from my data
%z=.025:.025:16; %change in depth in cm (or .25 mm)
subs=[];
for z=.025:.025:17
    index=round(z.*40);
    zf=z+.025;
    eo=(log(z/ap))/bp;
    ef=(log(zf/ap))/bp;
    s=((eo-ef)./(1+eo)).*(.025); %subsidence (cm) in each .025 mm area
    subs(index)=(s/1).*10000; %divide by 1 h (.025mm/h agg rate) to get subsidence rate (cm/h)*10000 to get microns per hour
end

depth=.025:.025:17;
subs=squeeze(subs);

%make subs cumulative
fs=flip(subs);
subs_cum=cumsum(fs);
subs_cum=flip(subs_cum);


%make the plot nice w/ subplots
h=figure
% Enlarge figure to full screen.
set(gcf,'Units', 'Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fontsize=6;

subplot(1,3,1)
Zmm=Z*10;
plot(e2,Zmm,'o','LineWidth',1)
axis ij
xlabel('void ratio (e)')
ylabel('depth (mm)')
title('Fig7a')
ylim([-30 170])
yticks([0 20 40 60 80 100 120 140 160])
%showfit exp
myfit1 = ezfit(e2,Zmm,'y(x) = a*exp(b*x)');
%myfit2 = ezfit(e3,Zmm,'y(x) = a*exp(b*x)');
%myfit3 = ezfit(e4,Zmm,'y(x) = a*exp(b*x)');
showfit(myfit1,'dispfitlegend','off','dispeqboxmode','off')
%showfit(myfit2,'dispfitlegend','off','dispeqboxmode','off')
%showfit(myfit3,'dispfitlegend','off','dispeqboxmode','off')
hold on

yline(0,'LineWidth',1, 'LineStyle', '--')

p=subplot(1,3,2)
hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick; 
%boxplot(subs2hr_mw,-10,'orientation','horizontal','Positions',-10,'Widths',10);
boxplot(subs10hr_mw,-10,'orientation','horizontal','Positions',-10,'Widths',10);
hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick;  
hold on
depthmm=depth*10;
%plot(subs_cumlow,depthmm,'LineWidth',2,'DisplayName','80% initial porosity')
plot(subs_cum,depthmm,'LineWidth',2,'DisplayName','90% initial porosity')
%plot(subs_cumhigh,depthmm,'LineWidth',2,'DisplayName','96% initial porosity') %add in curves for areas

%sometimes recieving marsh, sometimes not
xlabel('subsidence rate ({\mum}/h)')
%ylabel('depth (cm)')
ylim([-30 170])
xlim([0 350])
axis ij
set(gca,'Ytick',[0 20 40 60 80 100 120 140 160]);%force axis labels to be not messed up
set(gca,'YtickLabel',[0 20 40 60 80 100 120 140 160]);

legend('Location','southeast','FontSize',7);
yline(0,'LineWidth',1, 'LineStyle', '--','HandleVisibility','off')
text(20,-20,'Lidar measurements')

p=subplot(1,3,3)
hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick; 
%boxplot(subs2hr_mw,-10,'orientation','horizontal','Positions',-10,'Widths',10);
boxplot(subs10hr_mw,-10,'orientation','horizontal','Positions',-10,'Widths',10);
hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick;  
hold on
depthmm=depth*10;
plot(subs_cum,depthmm,'LineWidth',2,'DisplayName','100% marsh deposition')
plot(subs_cum.*.75,depthmm,'LineWidth',2,'DisplayName','75% marsh deposition')
plot(subs_cum.*.5,depthmm,'LineWidth',2,'DisplayName','50% marsh deposition') %add in curves for areas
plot(subs_cum.*.25,depthmm,'LineWidth',2,'DisplayName','25% marsh deposition')
%sometimes recieving marsh, sometimes not
xlabel('subsidence rate ({\mum}/h)')
%ylabel('depth (cm)')
ylim([-30 170])
xlim([0 350])
axis ij
set(gca,'Ytick',[0 20 40 60 80 100 120 140 160]);%force axis labels to be not messed up
set(gca,'YtickLabel',[0 20 40 60 80 100 120 140 160]);

legend('100% marsh deposition','75% marsh deposition','50% marsh deposition','25% marsh deposition','Location','southeast','FontSize',7);
yline(0,'LineWidth',1, 'LineStyle', '--','HandleVisibility','off')
text(20,-20,'Lidar measurements')
title('Fig7b')