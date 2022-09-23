%Temporal subsidence correlation-> how does subs at t1 correlate with subs
%at t2. CROPPED FOR ELEVATION

load ZD_19_2_dry
load subs19_2hr.mat
load marshmaps.mat
load subs19_10hr.mat


%separating each subs map into 10x10 squares for separate analysis through
%time. Square size was chosen bc it is generally smaller than max
%correlation length determined through variogram analysis

mean_cell_all19=[];
len_cell_all19=[];
std_cell_all19=[];
mean_cell_all19_marsh=[];
for i=1:280
    tstep=subs_cells_19_2hr(i);
    tstep_access=tstep{1};
    tstep_all_spots=tstep_access(:); %makes 1976x1 cell

    tstep_marsh=marsh_cells_19_2hr(i);
    tstep_access_marsh=tstep_marsh{1};
    tstep_all_spots_marsh=tstep_access_marsh(:); %makes 1976x1 cell
    
    for j=1:1976
        cell_1=cell2mat(tstep_all_spots(j));
        lencell=sum(~isnan(cell_1(:))); %number of non-nan values
        meancell=nanmean(cell_1(:));    %mean value of 10x10 chunk
        stdcell=nanstd(cell_1(:));
        %meancell(len_cell_1<40)==nan;
        if lencell<40               %if less than 40 pixels aren't nan, exclude
            meancell=nan;
        end
        %marsh values
        cell_marsh=cell2mat(tstep_all_spots_marsh(j));
        %lencell=sum(~isnan(cell_1(:))); %number of non-nan values
        meancell_marsh=nanmean(cell_marsh(:));    %mean value of 10x10 chunk
        %stdcell=nanstd(cell_1(:));
        
        len_cell_all19(j,i)=lencell;
        meancell(meancell==0)=nan;
        mean_cell_all19(j,i)=meancell;
        std_cell_all19_(j,i)=stdcell;
        
        mean_cell_all19_marsh(j,i)=meancell_marsh;
    end
end
%save cell_stats_19_zlim.mat mean_cell_all19_zlim len_cell_all19_zlim std_cell_all19_zlim mean_cell_all19_marsh_zlim
%for all timesteps
marsh_cells=mean_cell_all19_marsh(:);
subs_cells=mean_cell_all19(:);
figure
plot(marsh_cells,subs_cells,'o')
corrcoef(marsh_cells,subs_cells,'rows','complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_eachtime=[];
for q=1:280
    subscells_iter=mean_cell_all19(:,q);
    marshcells_iter=mean_cell_all19_marsh(:,q);
    
    %plot(marshcells_iter,subscells_iter,'o')
    R_eachtime(:,:,q)=corrcoef(marshcells_iter,subscells_iter,'rows','complete');
end
R_eachtime=squeeze(R_eachtime(2,1,:));

%do a simple linear regression for subs at t vs subs at t+2hr
R=[];
P=[];
for k=1:1976
    [R(:,:,k),P(:,:,k)] = corrcoef(mean_cell_all19(k,1:end-1),mean_cell_all19(k,2:end),'rows','complete');

end
R=squeeze(R(2,1,:));
R(R==1)=nan;
R(R==-1)=nan;
R19=R;
overall_R_19_2hr=nanmean(R19);
%%%%plotting stuff nicely%%%%
%%
h=figure
hold on
histogram(R_eachtime,'Normalization','probability','BinWidth',0.1) %marsh-subs
histogram(R19,'Normalization','probability','BinWidth',0.1) %subs-subs+2h corr
ylabel('probability')
xlabel('Pearson R')
xline(nanmedian(R19),'--','LineWidth',1)
xline(nanmedian(R_eachtime),'--','LineWidth',1)
legend({' 2 hour Marsh-{\sigma_s} Correlation','Correlation of {\sigma_s}(t) to {\sigma_s}(t+2h)'},'Location','northwest')
xlim([-1,1])
ylim([0,.27])
text(-.25,0.15,['Median = ',num2str(round(nanmedian(R19),2))],'FontSize',10)
text(0.47,0.25,['Median = ',num2str(round(nanmedian(R_eachtime),2))],'FontSize',10)
saveas(h,'marsh_subs_hists_FIG','pdf')





