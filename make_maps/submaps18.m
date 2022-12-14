%This script makes maps of subsidence at a each timestep for TDB-18.
%Ringing artifacts and areas receiving surface water flow are removed.
%Elevation screens can be placed on data be adjusting "lowbound_z" and
%"lowbound_z", then analyzing products with "..._zlim" at the end. 

%Outputs for further use: 
%    subs18_2hr.mat   -->  2 hour subsidence map
%    subs18_10hr.mat  -->  10 hour subsidence map

%load in data to workspace from figshare: 
load ZD.mat   %tdb18 dry Z data
load flowscreen18   %tdb18 flowscreens



%uncomment below to run filtering script
% dz_2hr_18_filt = zeros(796,522,280);
% for i = 1:2:561
%     index=(i+1)/2;
%     %index_hr=2.*index;
%     %each 2 hr increment
%     dz_2hr=ZD(:,:,i+2)-ZD(:,:,i);
%     %dz_2hr_map(:,:,index)=dz_2hr;
%     dz_2hr_18_filt(:,:,index)=Ringing_Script_3(dz_2hr);
% end
%save dz_2hr_18_filt.mat dz_2hr_18_filt

load dz_2hr_18_filt.mat

%initial parameters to set:
lowbound_z=0;  %low bound of elevation window in question
highbound_z=15; %high bound of elevation window in question
%max and min subsidence and marsh values to accept--> all values outside
%range will be set to nan:
max_subs_val=6000;
min_subs_val=-1000;


%stacking flowscreens into 2 hour screens (has flow been here at any point
%over 2 hours
%first timestep of ZD corresponds with hour1
subs18_2hr=zeros(796,522,280);
subs18_2hr_zlim=zeros(796,522,280);
flowscreen18_2hr=zeros(796,522,280);
for i = 1:2:559  
    index=(i+1)/2;
    flowscreen_2hr=flowscreen18(:,:,i).*flowscreen18(:,:,i+1).*flowscreen18(:,:,i+2);
    flowscreen18_2hr(:,:,index)=flowscreen_2hr;
    SL_end=25+0.25*(i+2);
    Z_min=SL_end+lowbound_z; %adding max and min elevations
    Z_max=SL_end+highbound_z;
    %Z_highbound=(ZD(:,:,i+2)<Z_max);
    %Z_lowbound=(ZD(:,:,i+2)>Z_min);   
    delta2=ZD(:,:,i+2);
    
    Z_highbound=(delta2<Z_max);
    Z_lowbound=(delta2>Z_min);
    RSLR_screen=delta2>SL_end;
    elev_screen=Z_highbound.*Z_lowbound; %elevation window
    zdiff_2hr_RSLR=dz_2hr_18_filt(:,:,index).*RSLR_screen; %screened for sea level at end of timestep
    
    %subs: 10 hour DOD created by summing together 10 1 hr DODs. Screened
    %for sea level at end of each 10 hour timestep and for surface water
    %flow at each hour.
    subs=-1000.*dz_2hr_18_filt(:,:,index).*flowscreen_2hr.*RSLR_screen;  %values in microns now, subsidence is + value, screened for flow
    subs(subs==0)=nan;
    subs(subs>max_subs_val)=nan;
    subs(subs<min_subs_val)=nan;
    subs18_2hr(:,:,index)=subs;%.*no_erosion;
    
    subs_zlim=-1000.*dz_2hr_18_filt(:,:,index).*flowscreen_2hr.*RSLR_screen.*elev_screen;  %values in microns now, subsidence is + value, screened for flow
    subs_zlim(subs_zlim==0)=nan;
    subs_zlim(subs_zlim>max_subs_val)=nan;
    subs_zlim(subs_zlim<min_subs_val)=nan;
    subs18_2hr_zlim(:,:,index)=subs_zlim;%.*no_erosion;  
    
end
save 'subs18_2hr.mat' 'subs18_2hr'
%save 'subs18_2hr_zlim.mat' 'subs18_2hr_zlim'

for j=2:2:280
   tstep4=j/2; 
   flowscreen_4hr=flowscreen18_2hr(:,:,j).*flowscreen18_2hr(:,:,j-1);
   subs18_4hr_zlim(:,:,tstep4)=(subs18_2hr_zlim(:,:,j)+subs18_2hr_zlim(:,:,j-1)).*flowscreen_4hr;
end

for k=3:3:280
   tstep6=k/3; 
   flowscreen_6hr=flowscreen18_2hr(:,:,k).*flowscreen18_2hr(:,:,k-1).*flowscreen18_2hr(:,:,k-2);
   subs18_6hr_zlim(:,:,tstep6)=(subs18_2hr_zlim(:,:,k)+subs18_2hr_zlim(:,:,k-1)+subs18_2hr_zlim(:,:,k-2)).*flowscreen_6hr;
end

for l=4:4:280
   tstep8=l/4; 
   flowscreen_8hr=flowscreen18_2hr(:,:,l).*flowscreen18_2hr(:,:,l-1).*flowscreen18_2hr(:,:,l-2).*flowscreen18_2hr(:,:,l-3);
   subs18_8hr_zlim(:,:,tstep8)=(subs18_2hr_zlim(:,:,l)+subs18_2hr_zlim(:,:,l-1)+subs18_2hr_zlim(:,:,l-2)+subs18_2hr_zlim(:,:,l-3)).*flowscreen_8hr;
end

for m=5:5:280
   tstep10=m/5; 
   flowscreen_10hr=flowscreen18_2hr(:,:,m).*flowscreen18_2hr(:,:,m-1).*flowscreen18_2hr(:,:,m-2).*flowscreen18_2hr(:,:,m-3).*flowscreen18_2hr(:,:,m-4);
   subs18_10hr_zlim(:,:,tstep10)=(subs18_2hr_zlim(:,:,m)+subs18_2hr_zlim(:,:,m-1)+subs18_2hr_zlim(:,:,m-2)+subs18_2hr_zlim(:,:,m-3)+subs18_2hr_zlim(:,:,m-4)).*flowscreen_10hr;
end

for n=5:5:280
   tstep10=n/5; 
   flowscreen_10hr=flowscreen18_2hr(:,:,n).*flowscreen18_2hr(:,:,n-1).*flowscreen18_2hr(:,:,n-2).*flowscreen18_2hr(:,:,n-3).*flowscreen18_2hr(:,:,n-4);
   subs18_10hr(:,:,tstep10)=(subs18_2hr(:,:,n)+subs18_2hr(:,:,n-1)+subs18_2hr(:,:,n-2)+subs18_2hr(:,:,n-3)+subs18_2hr(:,:,n-4)).*flowscreen_10hr;
end

% save subs18_4hr_zlim.mat subs18_4hr_zlim
% save subs18_6hr_zlim.mat subs18_6hr_zlim
% save subs18_8hr_zlim.mat subs18_8hr_zlim
% save subs18_10hr_zlim.mat subs18_10hr_zlim
save subs18_10hr.mat subs18_10hr

%make blocky subsidence maps by averaging over 10x10 pixel areas
subs18_2hr_blocky = zeros(796,522,280);
for i=1:280
    fun = @(block_struct)nanmean(nanmean((block_struct.data))) * ones(size(block_struct.data));
    submap=blockproc(subs18_2hr_zlim(:,:,i),[10 10],fun);
    subs18_2hr_blocky(:,:,i)=submap;
end
%save subs18_2hr_blocky.mat subs18_2hr_blocky

%%
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
SL180=25+0.25*180;
mm_dist=0:5:495;
h=figure
plot(mm_dist,prof170-SL180,'color','blue','LineWidth',2)
hold on
plot(mm_dist,prof180-SL180,'color','red','LineWidth',2)
xlabel('distance along transect T-T'' (mm)')
ylabel('elevation RSL at hour 180 (mm)')
pbaspect([2 1.4 1])
ylim([-7 12])
yline(0,'color','black','LineWidth',2)
text(230,8,'Hour 170','color','blue')
text(180,2.25,'Hour 180','color','red')
xl=[350 372.86];
yl=[10 8];
line(xl,yl,'color','black','LineWidth',2)
text(375,9.5,'5^{\circ} slope')
cd C:\Users\zapps\Box\Zapp\Manuscript\figs
saveas(h,'subs_scarp_FIG','pdf')