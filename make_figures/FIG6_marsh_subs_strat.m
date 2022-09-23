cd C:\Users\zapps\Box\Zapp\Subsidence_Maps_Procedure\input_matrices
load ZD_19_2_dry.mat
cd C:\Users\zapps\Box\Zapp\Subsidence_Maps_Procedure\output_maps
load marshmaps_zlim.mat
load subs19_2hr_zlim.mat
load subs18_2hr_zlim.mat
%load stratigraphy data
load mat_fraction.mat
load mat_thickness.mat

tot_marsh=nansum(marshmaps_zlim,3);
tot_marsh(tot_marsh==0)=NaN;

tot_subs=nansum(subs19_2hr_zlim,3);
tot_subs(tot_subs==0)=NaN;
tot_subs18=nansum(subs18_2hr_zlim,3);
tot_subs18(tot_subs18==0)=NaN;
%logical screen to determine when subsidence is being observed
subs_loc=zeros(750,747,280);
subs_loc18=zeros(796,522,280);
for i=1:280
    subs=subs19_2hr_zlim(:,:,i);
    subs=subs./2; %divide to make hourly
    subs_bool=~isnan(subs)==1;
    subs_loc(:,:,i)=subs_bool;
    subs18=subs18_2hr_zlim(:,:,i);
    subs_bool18=~isnan(subs18)==1;
    subs_loc18(:,:,i)=subs_bool18;   
end
subs_loc_tot=nansum(subs_loc,3);
subs_loc_tot(subs_loc_tot<50)=NaN; %remove areas in window less than 30 times
subs_avg=tot_subs./subs_loc_tot;

subs_loc_tot18=nansum(subs_loc18,3);
subs_loc_tot18(subs_loc_tot18<50)=NaN; %remove areas in window less than 30 times
subs_avg18=tot_subs18./subs_loc_tot18;
%subs and marsh as function of distance from back wall
idx=length(214:20:714);
subs_y=[];
marsh_y=[];
zthick_y=[];
for y=214:20:714
    totz19=ZD_19_2_dry(y:y+20,:,280); %modified 3/16/21 added : in 2nd dimension that was missing
    tm=tot_marsh(y:y+20,:);
    ts=subs_avg(y:y+20,:);
    zthick_y(y)=nanmean(nanmean(totz19));
    marsh_y(y)=nanmean(nanmean(tm));
    subs_y(y)=nanmean(nanmean(ts))./2;%divide by two to make hourly
end

subs_y=(subs_y)';
marsh_y=(marsh_y)';
zthick_y=(zthick_y)';

marsh_y_normed=(marsh_y./1000)./zthick_y;
%%%add in Jose's real data 2/22/21%%%%%
marshZ_downdip=nanmean(mat_thickness,2);
marshfrac_downdip=nanmean(mat_fraction,2);

ms=marshfrac_downdip(~isnan(marshfrac_downdip));
mf=marsh_y_normed(~isnan(marsh_y_normed));
msub=subs_y(subs_y~=0);
figure
plot(ms(1:17),mf(1:17),'o') %correlation between marsh fraction as function of dist from shoreline as calculated by maps & strat
corrcoef(ms(1:17),mf(1:17))
xlabel('marsh fraction (-) from stratigraphy')
ylabel('synthetic marsh fraction (-) from lidar maps')
hold on
xx=[0 0.6];
yy=[0 0.6];
line(xx,yy)
xlim([0 0.3])
ylim([0 0.6])
hold off
figure
plot(mf(1:17),msub(1:17),'o')
corrcoef(mf(1:17),msub(1:17))
xlabel('synthetic marsh fraction (-) from lidar maps')
ylabel('average subsidence (microns)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%see if i can interpolate marsh thicknesses

fill1=fillmissing(mat_thickness,'nearest',2);
fill2=fillmissing(fill1,'nearest',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
marshfrac_downdip=marshfrac_downdip(1:714);
subs_y(subs_y==0)=NaN;
marsh_y_normed(marsh_y_normed==0)=NaN;  %set zeros to nans
marsh_y_normed(553:end)=NaN;
subs_y(553:end)=NaN;

%subs_y, marsh_y_normed, marshfrac_downdip ->these are the arrays to use
pix=1:714;
mm_to_feed=(pix-214).*5;
idx_subs=~isnan(subs_y);
idx_marsh=~isnan(marsh_y_normed);
idx_marshstrat=~isnan(marshfrac_downdip);
%%
%%%start marsh-subs-strat fig
h=figure
% Enlarge figure to full screen.
set(gcf,'Units', 'Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fontsize=6;
subplot(3,3,[1,6])%%%%%%%%
plot(mm_to_feed(idx_marshstrat), marshfrac_downdip(idx_marshstrat),'s','color','#1f77b4','LineWidth',1)
hold on
ax = gca;
plot(mm_to_feed(idx_marsh),marsh_y_normed(idx_marsh),'-o','color','#1f77b4','LineWidth',1)  %mapped marsh fraction
ylabel('marsh fraction (-)','color','#1f77b4')
set(gca,'FontSize',fontsize)
xlabel('distance downstream (mm)')
set(gca,'FontSize',fontsize)
text(2400,0.025,'A')
yyaxis right
plot(mm_to_feed(idx_subs),subs_y(idx_subs),'-o','color','#FA6121','LineWidth',1)  %subsidence vs x direction
ylabel('average subsidence ({\mum}/h)','color','#FA6121')
set(gca,'FontSize',fontsize)
xline(1100,'LineWidth',2) %mean shoreline
ax.YAxis(1).Color = '#1f77b4';
ax.YAxis(2).Color = '#FA6121';
legend('stratigraphy','stacked marsh maps','subsidence rate', 'mean shoreline')
%text(mm_to_feed(idx_marsh),0.05,['A.'],'FontSize',8)

subplot(3,3,7)%%%%%%%%
plot(ms(1:17),mf(1:17),'o','MarkerSize',3,'color','#1f77b4') %correlation between marsh fraction as function of dist from shoreline as calculated by maps & strat
xlabel('marsh fraction (-) from stratigraphy')
set(gca,'FontSize',fontsize)
ylabel('marsh fraction (-) from lidar')
set(gca,'FontSize',fontsize)
hold on
Fit1=LinearModel.fit(ms(1:17),mf(1:17));
x1=0:0.05:0.6;
y1=(Fit1.Coefficients{2,1}.*x1)+(Fit1.Coefficients{1,1});
%axis equal
plot(x1,y1)
text(0.025,0.275,['R^2 = ',num2str(round(Fit1.Rsquared.Ordinary,2))],'FontSize',8)
text(0.25,0.035,'B')
% xx=[0 0.6];
% yy=[0 0.6];
% line(xx,yy,'color','black')
xlim([0 0.3])
ylim([0 0.3])
hold off
subplot(3,3,8)%%%%%%%%%
plot(mf(1:17),msub(1:17),'o','MarkerSize',3,'color','#1f77b4')
corrcoef(mf(1:17),msub(1:17))
xlabel('marsh fraction (-) from lidar')
set(gca,'FontSize',fontsize)
ylabel('average subsidence ({\mum})')
set(gca,'FontSize',fontsize)
xlim([0 0.3])
ylim([0 250])
Fit2=LinearModel.fit(mf(1:17),msub(1:17));
x2=0:0.05:0.6;
y2=(Fit2.Coefficients{2,1}.*x2)+(Fit2.Coefficients{1,1});
hold on
plot(x2,y2)
text(0.05,230,['R^2 = ',num2str(round(Fit2.Rsquared.Ordinary,3))],'FontSize',8)
text(0.25,30,'C')
hold off
subplot(3,3,9)%%%%%%
plot(ms(1:17),msub(1:17),'o','MarkerSize',3,'color','#1f77b4')
corrcoef(ms(1:17),msub(1:17))
xlabel('marsh fraction (-) from stratigraphy')
set(gca,'FontSize',fontsize)
ylabel('average subsidence ({\mum})')
set(gca,'FontSize',fontsize)
xlim([0 0.3])
ylim([0 250])
Fit3=LinearModel.fit(ms(1:17),msub(1:17));
x3=0:0.05:0.3;
y3=(Fit3.Coefficients{2,1}.*x3)+(Fit3.Coefficients{1,1});
hold on
plot(x3,y3)
text(0.025,230,['R^2 = ',num2str(round(Fit3.Rsquared.Ordinary,3))],'FontSize',8)
text(0.25,30,'D')
hold off
%LinearModel.fit(...)

saveas(h,'marsh_strat_subs_FIG','pdf')

%%









%%normalize marsh thickness by total thickness
%marsh_y_normed=marsh_y/zthick_y;

%same for control
% idx18=length(109:20:729);
% subs_y18=[];
% 
% for y=109:20:729
%     ts18=subs_avg18(y:y+20,:);
%     subs_y18(y)=nanmean(nanmean(ts18));
% end
% 
% subs_y18=(subs_y18)';


% sum each row of the matrix, then find rows with non-zero sum
idx_nonzerolines = sum(abs(subs_y),2)>0 ;
%idx_nonzerolines18 = sum(abs(subs_y18),2)>0 ;
% Create matrix B containing only the non-zero lines of A
subs_y = subs_y(idx_nonzerolines,:) ;
%subs_y18 = subs_y18(idx_nonzerolines18,:) ;
% sum each row of the matrix, then find rows with non-zero sum
idx_nonzerolines = sum(abs(marsh_y),2)>0 ;
% Create matrix B containing only the non-zero lines of A
marsh_y = marsh_y(idx_nonzerolines,:) ;
%marsh_y=marsh_y./10000; %micron to cm
zthick_y=zthick_y(idx_nonzerolines,:) ;

jose_idx=((linspace(1,26,26)).*100);
jose_idx_2=((linspace(1,17,17)).*100);

marsh_y_normed=marsh_y./(zthick_y*1000); %marsh fraction from stacked marsh maps
%%%%
%make jose's data
dist_downstream= [1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2400 2500 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700];
marsh_fraction= [.03 .05 .075 .07 .08 .115 .135 .175 .22 .23 .22 .19 .135 .13 .128 .12 .128 .129 .145 .135 .08 .085 .04];
%plot marsh frac from strat, marsh frac from maps, avg subs as function if
%dist basinward
figure
plot(dist_downstream, marshfrac_downdip,'s','color','blue')
hold on
%yyaxis left
ax = gca;
plot(jose_idx(1:17),marsh_y_normed(1:17),'-o','color','blue','LineWidth',2)  %mapped marsh fraction
ylabel('marsh fraction (-)','color','blue')
xlabel('distance downstream (mm)')


yyaxis right
plot(jose_idx_2,subs_y,'-o','color','#FA6121','LineWidth',2)  %subsidence vs x direction
ylabel('average subsidence (microns/h)','color','#FA6121')

xline(2200,'LineWidth',2) %mean shoreline
ax.YAxis(1).Color = 'blue';
ax.YAxis(2).Color = '#FA6121';
legend('strat','marshmaps','subsidence', 'mean shoreline')
%add in control subs trend ???
%plot(jose_idx_2(1:end-1),subs_y18,'--','color','b')

%set(gcf,'color','w');
%cd C:\Users\zapps\Downloads\addaxis6
%figure
%plot(dist_downstream, marsh_fraction,'s','color','magenta')
%addaxis(jose_idx,marsh_y)
%addaxis(jose_idx_2,subs_y)
%addaxislabel(1,'one');
%addaxislabel(2,'two');
%addaxislabel(3,'three');
%xlabel('distance downstream (mm)')
%ylabel('marsh fraction (-)')
%jose_idx=((linspace(1,26,26)).*100)+1000;
%jose_idx_2=((linspace(1,17,17)).*100)+1000;


%legend('marsh thickness (micron)',('avg subs (micron)')
%std to determine similarity
devsubs=nanstd(subs_avg(:));
devmarsh=nanstd(tot_marsh(:));
subs_avg_dev=subs_avg./devsubs;  %number of standard deviations for each point
tot_marsh_dev=tot_marsh./devmarsh;
marsh_subs_std_sim=1./(subs_avg_dev./tot_marsh_dev);
imagesc(marsh_subs_std_sim)
caxis([0,2.5])
colorbar
figure
%mean to determine similarity
meanstatsubs=nanmean(subs_avg(:));
meanstatmarsh=nanmean(tot_marsh(:));
subs_avg_meanstat=subs_avg./meanstatsubs;  %value of each subs pixel normalized by mean
tot_marsh_meanstat=tot_marsh./meanstatmarsh;  %value of each marsh pixel normalized by mean
marsh_subs_mean_sim=1./(subs_avg_meanstat./tot_marsh_meanstat);   %similarity score
imagesc(marsh_subs_mean_sim)
caxis([0,2.5])
colorbar
%why not use zscore you dummy
zscoresubs=(subs_avg-nanmean(subs_avg(:)))./nanstd(subs_avg(:));
zscoremarsh=(tot_marsh-nanmean(tot_marsh(:)))./nanstd(tot_marsh(:));
%subs_avg_dev=subs_avg./devsubs;  %number of standard deviations for each point
%tot_marsh_dev=tot_marsh./devmarsh;
marsh_subs_zscore_sim=1./(abs(zscoresubs-zscoremarsh));
imagesc(marsh_subs_zscore_sim)
caxis([0,2.5])
colorbar

%%%
plot(jose_idx_2(1:end-1),subs_y18,'--')
