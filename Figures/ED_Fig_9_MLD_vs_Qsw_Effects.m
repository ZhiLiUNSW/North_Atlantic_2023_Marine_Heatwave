%% Extended Data Figure: Map of Dominant Terms for MLT Bbudget Anomaly in 2023 Relative to 1981-2010
%  June, July 
%  Monthly Fields from IAP and ERA5, 1981-2023


%% ########################################################################
%%% Plotting 1: Dominant Terms for MLT Bbudget Anomaly in 2023 Relative to 1981-2010
clc;clear
time_ann=(1981:2022)';


% #########################################################################
% #########################################################################
% 1. MLTa from IAP data
    % Monthly MLT during 1981-2023
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA(:,:,:,count)=mlt_monthly(:,:,1:12);
          clear mlt_monthly
          count=count+1;
    end
    clear year count Sxy_NA lon lat

    % MLT tendency
    dmlt_monthly_NA = mlt_monthly_NA(:,:,2:12,:) - mlt_monthly_NA(:,:,1:11,:); % mon: 1.5:11.5
    % Anomaly relative to the clim mean
    dmlt_monthly_2023 = dmlt_monthly_NA(:,:,:,43) - squeeze(nanmean(dmlt_monthly_NA(:,:,:,1:30),4));
    % Centred at mon: 2:11
    dmlt_monthly_2023_C(1:size(dmlt_monthly_2023,1),1:size(dmlt_monthly_2023,2),1)    = NaN;
    dmlt_monthly_2023_C(1:size(dmlt_monthly_2023,1),1:size(dmlt_monthly_2023,2),2:11) = (dmlt_monthly_2023(:,:,1:10)+dmlt_monthly_2023(:,:,2:11))./2;
    dmlt_monthly_2023_C(1:size(dmlt_monthly_2023,1),1:size(dmlt_monthly_2023,2),  12) = NaN;
    clear mlt_month*  dmlt_monthly_2023 dmlt_monthly_NA
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2. MLT Budget in 1981-2023 - MLD' Term
count=1;
for year=1981:2023
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'_2_MLD_prime.mat'],...
          'MLT_Qswr_mon')
    dmlt_Qswr_mon0(:,:,:,count)=MLT_Qswr_mon(:,:,1:12);
    clear MLT_Qswr_mon
    count=count+1;
end
    % Anomaly relative to the clim mean
    dmlt_Qswr_2023_1_MLD_prime = dmlt_Qswr_mon0(:,:,:,43) - squeeze(nanmean(dmlt_Qswr_mon0(:,:,:,1:30),4));
    clear dmlt_Qswr_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 3. MLT Budget in 1981-2023 - Qsw' Term
count=1;
for year=1981:2023
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'_3_Qsw_prime.mat'],...
          'MLT_Qswr_mon')
    dmlt_Qswr_mon0(:,:,:,count)=MLT_Qswr_mon(:,:,1:12);
    clear MLT_Qswr_mon
    count=count+1;
end
    % Anomaly relative to the clim mean
    dmlt_Qswr_2023_2_Qswr_prime = dmlt_Qswr_mon0(:,:,:,43) - squeeze(nanmean(dmlt_Qswr_mon0(:,:,:,1:30),4));
    clear dmlt_Qswr_mon0
% #########################################################################
% #########################################################################

% 4. Residual terms:

% 5. Dominant terms
dmlt_2023_domin_ratio = dmlt_Qswr_2023_2_Qswr_prime - dmlt_Qswr_2023_1_MLD_prime;
% #########################################################################
% #########################################################################



%% ########################################################################
%  Figure 2.3: June, July, Difference between MLD and QSW term
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.080; ixe = 0.150;  ixd = 0.04; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.220; iye = 0.200;  iyd = 0.06; iyw = (1-iys-iye-0*iyd)/1;

pos{11}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:5,:)  = (color(13:-1:9,:).*2/3+color(14:-1:10,:).*1/3)./1;
color0(6:10,:) = (color( 7:-1:3,:).*1/3+color( 6:-1: 2,:).*2/3)./1;
% #########################################################################


% #########################################################################
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(dmlt_2023_domin_ratio(:,:,6)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',23,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-2 2]);
        title('a. Qsw vs MLD effects (June 2023)','fontsize',24,'FontWeight','bold')
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(dmlt_2023_domin_ratio(:,:,7)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',23,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-2 2]);
        title('b. Qsw vs MLD effects (July 2023)','fontsize',24,'FontWeight','bold')
        
        

        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.016 iys+0.1*iyw+0*iyd 0.012 0.8*iyw+0*iyd]);
        set(hBar1, 'ytick',-2:0.4:2,'yticklabel',{'<-2','-1.6','-1.2','-0.8','-0.4','0','0.4','0.8','1.2','1.6','>2'},...
                   'fontsize',23,'FontName','Arial','LineWidth',1.2,'TickLength',0.045);
        
        ylabel(hBar1, '[ \circC per month ]','rotation',90); 
    
% #########################################################################



% #########################################################################
% #########################################################################
   
disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')



