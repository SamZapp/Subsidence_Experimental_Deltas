%subs profile comparison
load ZD.mat
load ZD_19_2_dry.mat
load subs19_10hr.mat
%get transect for scarp analysis hr 170-180 (3/18/21)
topo160=ZD(:,:,160);
topo170=ZD(:,:,170);
topo180=ZD(:,:,180);
topo190=ZD(:,:,190);

xi=[404 404];
yi=[355 454];
prof160=improfile(topo160,xi,yi);
prof170=improfile(topo170,xi,yi);
prof180=improfile(topo180,xi,yi);
prof190=improfile(topo190,xi,yi);
SL160=25+0.25*160;
SL190=25+0.25*190;
mm_dist=0:5:495;

h=figure
%h=subplot(1,2,1)
h1=plot(mm_dist,prof160-SL160,'LineWidth',2)
hold on
h2=plot(mm_dist,prof170-SL160,'LineWidth',2)
h3=plot(mm_dist,prof180-SL160,'LineWidth',2)
h4=plot(mm_dist,prof190-SL160,'LineWidth',2)
% plot(mm_dist,prof170-SL180,'color','blue','LineWidth',2)
% hold on
% plot(mm_dist,prof180-SL180,'color','red','LineWidth',2)
%legend([h1 h2 h3 h4],'hour 160','hour 170','hour 180','hour 190')
xlabel('distance along transect C-C'' (mm)')
ylabel('elevation RSL at hour 160 (mm)')
pbaspect([2 1.4 1])
ylim([-7 20])
yline(0,'color','black','LineWidth',2)
yline(.25*30,'--','color','black','LineWidth',2)
% text(230,8,'Hour 170','color','blue')
% text(180,2.25,'Hour 180','color','red')
xl=[350 372.86];
yl=[10 8];
line(xl,yl,'color','black','LineWidth',2)
text(375,9.5,'5^{\circ} slope')
legend([h1 h2 h3 h4],'hour 160','hour 170','hour 180','hour 190')
text(25,-1,'RSL at hr 160')
text(25,6.5,'RSL at hr 190')
hold off
%%%treatment
topo500=ZD_19_2_dry(:,:,250);
topo510=ZD_19_2_dry(:,:,255);
topo520=ZD_19_2_dry(:,:,260);
topo530=ZD_19_2_dry(:,:,265);

xi=[454 454];
yi=[395+80 494+80];
prof500=improfile(topo500,xi,yi);
prof510=improfile(topo510,xi,yi);
prof520=improfile(topo520,xi,yi);
prof530=improfile(topo530,xi,yi);
SL500=25+0.25*500;
SL530=25+0.25*530;
mm_dist=0:5:495;
%h=subplot(1,2,2)
h=figure
h5=plot(mm_dist,prof500-SL500,'LineWidth',2)
hold on

h6=plot(mm_dist,prof510-SL500,'LineWidth',2)
h7=plot(mm_dist,prof520-SL500,'LineWidth',2)
h8=plot(mm_dist,prof530-SL500,'LineWidth',2)
% plot(mm_dist,prof170-SL180,'color','blue','LineWidth',2)
% hold on
% plot(mm_dist,prof180-SL180,'color','red','LineWidth',2)
xlabel('distance along transect T-T'' (mm)')
ylabel('elevation RSL at hour 500 (mm)')
pbaspect([2 1.4 1])
ylim([-7 20])
yline(0,'color','black','LineWidth',2)
yline(.25*30,'--','color','black','LineWidth',2)
% text(230,8,'Hour 170','color','blue')
% text(180,2.25,'Hour 180','color','red')
xl=[350 372.86];
yl=[10 8];
line(xl,yl,'color','black','LineWidth',2)
text(375,9.5,'5^{\circ} slope')
legend([h5 h6 h7 h8],'hour 500','hour 510','hour 520','hour 530')
text(25,-1,'RSL at hr 500')
text(200,8.5,'RSL at hr 530')     %455
quiver(140,8.5,0,-subs19_10hr(503,454,53)/1000,'LineWidth',3,'color','#7E2F8E')
quiver(260,6.5,0,-subs19_10hr(527,454,52)/1000,'LineWidth',3,'color','#EDB120')
quiver(325,3.8,0,-subs19_10hr(540,454,51)/1000,'LineWidth',3,'color','#D95319')
% cd C:\Users\zapps\Box\Zapp\Manuscript\figs
% saveas(h,'subs_scarp_FIG','pdf')
subsprof510=improfile(subs19_10hr(:,:,51),xi,yi);
subsprof520=improfile(subs19_10hr(:,:,52),xi,yi);
subsprof530=improfile(subs19_10hr(:,:,53),xi,yi);
figure
plot(mm_dist,subsprof510,'LineWidth',2)