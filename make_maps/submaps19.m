%This script makes maps of subsidence and marsh deposit thickness at a 2hr
%timestep for TDWB-19-2. These are made using "dry" (beginning of every
%other run hour), and "wet" (48 minutes into each run hour) lidar scans. 
%Step 1: make difference maps of subsidence at each 2 hour interval (dry(t+1)-dry(t))
%and initial marsh deposit thickness (wet(t)-dry(t)). These difference maps 
%are passed through a filter (Ringing_sScript_3) which mitigates ringing artifacts.
%Step2: We then filter the difference maps to remove areas outside of the desired 
%elevation range,areas recieving fluvial reworking (ie areas with surface water),
%and areas with outlier values. 
%Elevation screens can be placed on data be adjusting "lowbound_z" and
%"lowbound_z", then analyzing products with "..._zlim" at the end. 

%Outputs for further use: 
%    subs19_2hr.mat   -->  2 hour subsidence map
%    marshmaps.mat    -->  initial marsh deposit map
%    subs19_10hr_mw.mat    -->  10 hour subsidence within marsh window for comparison in Fig7b   


%load in data to workspace from figshare: 
load 'ZD_19_2_dry.mat'  %topography at every even hour of runtime
load 'ZD_19_2_wet.mat'  %topography 48 minutes into each hour  
load 'flowscreen19.mat'  %screen for surface water based on color threshold on wet scans


%uncomment to go through procedure for mitigating ringing artifacts
% dz_2hr_filt = zeros(750,747,280);
% dz_marsh_filt = zeros(750,747,280);
% for i = 1:2:559
%     index=(i+1)/2;
%     %maps of marsh deposits
%     dz_marsh=ZD_19_2_wet(:,:,i)-ZD_19_2_dry(:,:,index);
%     %each 2 hr increment
%     dz_2hr=(ZD_19_2_dry(:,:,index+1)-ZD_19_2_wet(:,:,i));
%     %dz_2hr_map(:,:,index)=dz_2hr;
%     dz_marsh_filt(:,:,index)=Ringing_Script_3(dz_marsh);
%     dz_2hr_filt(:,:,index)=Ringing_Script_3(dz_2hr);
% end
% save TDWB_19_2_filtered_dz_maps.mat dz_marsh_filt dz_2hr_filt

load TDWB_19_2_filtered_dz_maps.mat

%initial parameters to set:
lowbound_z=0;  %low bound of elevation window in question
highbound_z=15; %high bound of elevation window in question
highbound_z_marsh=5; %highbound for marsh map-> this is upper bound of marsh window
%max and min subsidence and marsh values to accept--> all values outside
%range will be set to nan:
max_subs_val=6000;
min_subs_val=-1000;
max_marsh_val=6000;
min_marsh_val=-1000;

%make 2 hr subs maps
subs19_2hr_zlim = zeros(750,747,280); %preallocate subsidence maps
marshmaps_zlim = zeros(750,747,280);  %preallocate marsh deposit maps
subs19_2hr = zeros(750,747,280);
marshmaps = zeros(750,747,280);
subs19_2hr_activemarsh = zeros(750,747,280);
flowscreen19_2hr = zeros(750,747,280);
for i = 1:2:559
    index=(i+1)/2;
    %screen out surface water by color threshold
    flowscreen_2hr=flowscreen19(:,:,i).*flowscreen19(:,:,i+1); % screens out flow based on color threshold of both wet scans in timestep
    flowscreen19_2hr(:,:,index)=flowscreen_2hr;
    marshmap_flowscreen=flowscreen19(:,:,i); %screens out flow for 
    
    SL_end=25+0.25*(i+1); %SL at hour 2,4,6...
    Z_min=SL_end+lowbound_z; %adding max and min elevations
    Z_max=SL_end+highbound_z;
    Z_max_marsh=SL_end+highbound_z_marsh;
    dry_t=i+1; %2,3,4,5->hr 2,4,6,8  (bc we start at t=0 here)
    Z_highbound=(ZD_19_2_dry(:,:,index+1)<Z_max);
    Z_lowbound=(ZD_19_2_dry(:,:,index+1)>Z_min);
    Z_highbound_marsh=(ZD_19_2_dry(:,:,index+1)<Z_max_marsh);
    elev_screen=Z_highbound.*Z_lowbound; %elevation window
    elev_screen_marsh=Z_highbound_marsh.*Z_lowbound; %elevation window
    %zscreen(:,:,index)=RSLR_screen;
    RSLR_screen=ZD_19_2_dry(:,:,index+1)>SL_end;
    %make sea level band for plotting
    %SL_band_high=SL_end+0.5;
    %SL_band_low=SL_end-0.5;
    SL_band_19(:,:,index)=ZD_19_2_dry(:,:,index+1)>SL_end-0.25 & ZD_19_2_dry(:,:,index+1)<SL_end+0.25;
    
    subs=-1000.*flowscreen_2hr.*RSLR_screen.*(10/6).*(dz_2hr_filt(:,:,index)); %subs in microns (and is positive). !!!note 10/6 correction for time
    subs(subs==0)=nan;
    subs(subs>max_subs_val)=nan;
    subs(subs<min_subs_val)=nan;
    index=(i+1)/2;
    subs19_2hr(:,:,index)=subs;
    marsh=1000.*marshmap_flowscreen.*RSLR_screen.*(dz_marsh_filt(:,:,index));  %didnt screen for SL for marsh maps
    marsh(marsh==0)=nan;
    marsh(marsh>max_marsh_val)=nan; %this is roughly the largest value for any of the sediment mounds we accidentally created
    marsh(marsh<min_marsh_val)=nan;  %measurement error
    marshmaps(:,:,index)=marsh;
    
    
    subs_zlim=-1000.*flowscreen_2hr.*elev_screen.*(10/6).*(dz_2hr_filt(:,:,index)); %subs in microns (and is positive). !!!note 10/6 (or 120/72) correction for lost time
    subs_zlim(subs_zlim==0)=nan;
    subs_zlim(subs_zlim>max_subs_val)=nan;
    subs_zlim(subs_zlim<min_subs_val)=nan;
    index=(i+1)/2;
    subs19_2hr_zlim(:,:,index)=subs_zlim;
    
    subs_active=-1000.*flowscreen_2hr.*elev_screen.*(10/6).*(dz_2hr_filt(:,:,index));  %didnt screen for SL for marsh maps
    subs_active(subs_active==0)=nan;
    subs_active(subs_active>max_marsh_val)=nan; %this is roughly the largest value for any of the sediment mounds we accidentally created
    subs_active(subs_active<min_marsh_val)=nan;  %measurement error
    subs19_2hr_activemarsh(:,:,index)=subs_active;
    
    marshs=1000.*marshmap_flowscreen.*elev_screen_marsh.*(dz_marsh_filt(:,:,index));  %didnt screen for SL for marsh maps
    marshs(marshs==0)=nan;
    marshs(marshs>max_marsh_val)=nan; %this is roughly the largest value for any of the sediment mounds we accidentally created
    marshs(marshs<min_marsh_val)=nan;  %measurement error
    marshmaps_zlim(:,:,index)=marshs;
end


%save subs19_2hr_zlim.mat subs19_2hr_zlim -v7.3
save marshmaps_zlim.mat marshmaps_zlim
save subs19_2hr.mat subs19_2hr -v7.3
save marshmaps.mat marshmaps
%save SL_band_19.mat SL_band_19


for j=2:2:280
   tstep4=j/2;
   flowscreen_4hr=flowscreen19_2hr(:,:,j).*flowscreen19_2hr(:,:,j-1);
   subs19_4hr_zlim(:,:,tstep4)=(subs19_2hr_zlim(:,:,j)+subs19_2hr_zlim(:,:,j-1)).*flowscreen_4hr;
end

for k=3:3:280
   tstep6=k/3;
   flowscreen_6hr=flowscreen19_2hr(:,:,k).*flowscreen19_2hr(:,:,k-1).*flowscreen19_2hr(:,:,k-2);
   subs19_6hr_zlim(:,:,tstep6)=(subs19_2hr_zlim(:,:,k)+subs19_2hr_zlim(:,:,k-1)+subs19_2hr_zlim(:,:,k-2)).*flowscreen_6hr;
end

for l=4:4:280
   tstep8=l/4;
   flowscreen_8hr=flowscreen19_2hr(:,:,l).*flowscreen19_2hr(:,:,l-1).*flowscreen19_2hr(:,:,l-2).*flowscreen19_2hr(:,:,l-3);
   subs19_8hr_zlim(:,:,tstep8)=(subs19_2hr_zlim(:,:,l)+subs19_2hr_zlim(:,:,l-1)+subs19_2hr_zlim(:,:,l-2)+subs19_2hr_zlim(:,:,l-3)).*flowscreen_8hr;
end

for m=5:5:280
   tstep10=m/5;
   flowscreen_10hr=flowscreen19_2hr(:,:,m).*flowscreen19_2hr(:,:,m-1).*flowscreen19_2hr(:,:,m-2).*flowscreen19_2hr(:,:,m-3).*flowscreen19_2hr(:,:,m-4);
   subs19_10hr_zlim(:,:,tstep10)=(subs19_2hr_zlim(:,:,m)+subs19_2hr_zlim(:,:,m-1)+subs19_2hr_zlim(:,:,m-2)+subs19_2hr_zlim(:,:,m-3)+subs19_2hr_zlim(:,:,m-4)).*flowscreen_10hr;
end

for n=5:5:280
   tstep10=n/5; 
   flowscreen_10hr=flowscreen19_2hr(:,:,n).*flowscreen19_2hr(:,:,n-1).*flowscreen19_2hr(:,:,n-2).*flowscreen19_2hr(:,:,n-3).*flowscreen19_2hr(:,:,n-4);
   subs19_10hr(:,:,tstep10)=(subs19_2hr(:,:,n)+subs19_2hr(:,:,n-1)+subs19_2hr(:,:,n-2)+subs19_2hr(:,:,n-3)+subs19_2hr(:,:,n-4)).*flowscreen_10hr;
end

%subsidence in marsh window for FIG7b
for i=1:280
    s=subs19_2hr_activemarsh(:,:,i);
    subs2hr_mw(i)=nanmean(nanmean(s));
end
s10=reshape(subs2hr_mw,5,56);
subs19_10hr_mw=sum(s10)/10;
save subs19_10hr_mw.mat subs19_10hr_mw

% save subs19_4hr_zlim.mat subs19_4hr_zlim
% save subs19_6hr_zlim.mat subs19_6hr_zlim
% save subs19_8hr_zlim.mat subs19_8hr_zlim
% save subs19_10hr_zlim.mat subs19_10hr_zlim
save subs19_10hr.mat subs19_10hr

%make blocky subsidence maps by averaging over 10x10 pixel areas
subs19_2hr_blocky = zeros(750,747,280);
marshmaps_blocky = zeros(750,747,280);
for i=1:280
    fun = @(block_struct)nanmean(nanmean((block_struct.data))) * ones(size(block_struct.data));
    submap=blockproc(subs19_2hr_zlim(:,:,i),[10 10],fun);
    marshmap=blockproc(marshmaps_zlim(:,:,i),[10 10],fun);
    subs19_2hr_blocky(:,:,i)=submap;
    marshmaps_blocky(:,:,i)=marshmap;
end
save subs19_2hr_blocky.mat subs19_2hr_blocky
save marshmaps_blocky.mat marshmaps_blocky




