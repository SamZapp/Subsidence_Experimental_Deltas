%reproduce Fig 3d

%load in data to workspace from figshare: 
load subs19_2hr.mat
load subs18_2hr.mat

load subs19_10hr.mat
load subs18_10hr.mat

load ZD.mat
load ZD_19_2_dry.mat

scarpZmean18_tot=zeros(280,50);
scarpZmean19_tot=zeros(280,50);
for i=1:280
    RSL18=25+i*0.5;
    RSL19=25+i*0.5;
    %RSL=25+i*0.5;
    ZD18=ZD(:,:,2*i-1); ZRSL18=ZD18-RSL18;
    ZD19=ZD_19_2_dry(:,:,i+1); ZRSL19=ZD19-RSL19;
    subs18=subs18_2hr(:,:,i);
    scarpZmean18=zeros(1,50);
    subs19=subs19_2hr(:,:,i);
    scarpZmean19=zeros(1,50);
    for z=1:50
        upbound18=(ZRSL18<=z);
        lowbound18=(ZRSL18>z-1);
        scarp18=upbound18.*lowbound18;
        upbound19=(ZRSL19<=z);
        lowbound19=(ZRSL19>z-1);
        scarp19=upbound19.*lowbound19;
        
        scarpZ18=subs18.*scarp18;
        scarpZ18(scarpZ18==0)=NaN;
        scarpZmean18(z)=nanmean(scarpZ18,'all');
        scarpZ19=subs19.*scarp19;
        scarpZ19(scarpZ19==0)=NaN;
        scarpZmean19(z)=nanmean(scarpZ19,'all');        
    end
    scarpZmean18_tot(i,:)=scarpZmean18;
    scarpZmean19_tot(i,:)=scarpZmean19;
end

subs_profile18=nanmean(scarpZmean18_tot)/2;
subs_profile19=nanmean(scarpZmean19_tot)/2;

figure
%set(gca,'Units','inches','Position',[0.5 0.5 2.25 2],'FontSize',12)
plot(1:50,subs_profile18,'LineWidth',2)
hold on
plot(1:50,subs_profile19,'LineWidth',2)
xlim([0 25])
ylim([0 350])
xlabel('elevation above sea level (mm)','FontSize',12)
ylabel('subsidence rate ({\mum}/hr)','FontSize',12)
legend('Control','Treatment')
rectangle('Position',[0,0,5,350],'FaceColor',[0.8500 0.3250 0.0980 0.3],'EdgeColor','none')
