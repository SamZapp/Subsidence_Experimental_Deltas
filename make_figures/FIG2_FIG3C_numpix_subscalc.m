%Reproduce Fig2 and Fig3c


%load in data to workspace from figshare: 
load subs19_2hr.mat
load subs18_2hr.mat
load subs19_10hr.mat
load subs18_10hr.mat

npix18=[];
npix19=[];
for i=1:280
    s19=subs19_2hr(:,:,i);
    npix19(i)=length(s19(~isnan(s19)));
    s18=subs18_2hr(:,:,i);
    npix18(i)=length(s18(~isnan(s18)));
end
szn19=1000.*(sqrt(((1./npix19).*sqrt(1)).^2)); %uncertainty as a function of number of pixels



npix18_10hr=[];
npix19_10hr=[];
for i=1:56
    s19=subs19_10hr(:,:,i);
    npix19_10hr(i)=length(s19(~isnan(s19)));
    s18=subs18_10hr(:,:,i);
    npix18_10hr(i)=length(s18(~isnan(s18)));
end

hr=2:2:560;
hr10=10:10:560;
dx=5/1000; %m ->pixel width

%average subsidence
%load
s18avg=[];
s19avg=[];
s18avg_10hr=[];
s19avg_10hr=[];
for i=1:56
    s18_10=subs18_10hr(:,:,i);
    s19_10=subs19_10hr(:,:,i);
%     s18_2=subs18_2hr(:,:,i);
%     s19_2=subs19_2hr(:,:,i);
    s18c=s18_10(s18_10>-1000 & s18_10<5000);
    s19c=s19_10(s19_10>-1000 & s19_10<5000);
%     s18c2=s18_2(s18_2>-1000 & s18_2<5000);
%     s19c2=s19_2(s19_2>-1000 & s19_2<5000);
    s18c=s18c/10;
    s19c=s19c/10;
%     s18c2=s18c2/2;
%     s19c2=s19c2/2;
    s18avg(i)=nanmean(nanmean(s18c));
    s19avg(i)=nanmean(nanmean(s19c));
%     s18avg2hr(i)=nanmean(nanmean(s18c2));
%     s19avg2hr(i)=nanmean(nanmean(s19c2));
end  

for i=1:280
%     s18_10=subs18_10hr(:,:,i);
%     s19_10=subs19_10hr(:,:,i);
    s18_2=subs18_2hr(:,:,i);
    s19_2=subs19_2hr(:,:,i);
%     s18c=s18_10(s18_10>-1000 & s18_10<5000);
%     s19c=s19_10(s19_10>-1000 & s19_10<5000);
    s18c2=s18_2(s18_2>-1000 & s18_2<5000);
    s19c2=s19_2(s19_2>-1000 & s19_2<5000);
%     s18c=s18c/10;
%     s19c=s19c/10;
    s18c2=s18c2/2;
    s19c2=s19c2/2;
%     s18avg(i)=nanmean(nanmean(s18c));
%     s19avg(i)=nanmean(nanmean(s19c));
    s18avg2hr(i)=nanmean(nanmean(s18c2));
    s19avg2hr(i)=nanmean(nanmean(s19c2));
end 
addnan=(ones(280-56,1)*NaN)';
s18_avg=[s18avg addnan];
s19_avg=[s19avg addnan];
subsavg=[s18avg2hr; s18_avg; s19avg2hr; s19_avg]';

%make FIG2
figure
subplot(1,2,1)
set(gca,'Units','inches','Position',[0.5 0.5 2.25 2],'FontSize',12)
plot(hr,npix18*(dx^2),'.','LineWidth',1,'color','#0072BD')
axis square
hold on
plot(hr,npix19*(dx^2),'.','LineWidth',1,'color','#D95319')
plot(hr10,npix18_10hr*(dx^2),'LineWidth',2,'color','#0072BD')
plot(hr10,npix19_10hr*(dx^2)','LineWidth',2,'color','#D95319')
xlabel('run time (hr)','FontSize',12)
%ylabel('number of measurements')
ylabel('Area Measured (m^2)','FontSize',12)
%yline(10000*(dx^2),'-','LineWidth',2);%10 micron thresh
yline(1600*(dx^2),'-','LineWidth',2);%25 micron thresh
legend('Control (2 hr)','Treatment (2hr)','Control (10 hr)','Treatment (10 hr)','25 micron error threshold')
ylim([0 3])
set(gca,'Ydir','normal')

subplot(1,2,2)
set(gca,'Units','inches','Position',[3 0.5 2.25 2],'FontSize',12)
boxplot(subsavg,'Labels',{'2 hr','10 hr','2 hr','10 hr'})
ylabel('Average Subsidence Rate ({\mum}/hr)','FontSize',12) 
yline(250,'-','LineWidth',2);
axis square
set(gca,'Ydir','normal')

%make FIG3C
figure                
plot(10:10:560,s18avg,'LineWidth',2)
hold on
plot(10:10:560,s19avg,'LineWidth',2)
xlabel('run time (hr)','FontSize',12)
ylabel('subsidence rate (\mu/hr)','FontSize',12)
