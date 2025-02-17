%% Extended Data Figure: Monthly Cloud Cover Anomaly in 2023 from 1981-2010 Mean
%  Low and Medium Cloud Cover
%  May, June, July, August, Monthly Mean 1981-2023
%  ERA-5 Monthly Averaged Reanalysis on Single Levels


%% #######################################################################
%% Figure #1: 2023 Cloud Cover Anomaly Relative to 1981-2010 Mean
clc;clear

% #########################################################################
% #########################################################################
% 1. High, Medium, and Low Cloud Cover Climatology during 1981-2010
    count=1;
    for year = 1981:2010
          disp(['ERA5 Monthly HML Cloud Year# ',num2str(year)])
          load(['cloud_hml_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'hcloud_total','mcloud_total','lcloud_total','lon','lat')
          hcloud_mon0(:,:,:,count)=squeeze(hcloud_total(:,1:end,:));
          mcloud_mon0(:,:,:,count)=squeeze(mcloud_total(:,1:end,:));
          lcloud_mon0(:,:,:,count)=squeeze(lcloud_total(:,1:end,:));
          clear hcloud_total mcloud_total lcloud_total
          count=count+1;
    end
    clear year count
    % high cloud
    hcloud_mon_clim_May=squeeze(nanmean(nanmean(hcloud_mon0(:,:,5,:),4),3));
    hcloud_mon_clim_Jun=squeeze(nanmean(nanmean(hcloud_mon0(:,:,6,:),4),3));
    hcloud_mon_clim_Jul=squeeze(nanmean(nanmean(hcloud_mon0(:,:,7,:),4),3));
    hcloud_mon_clim_Aug=squeeze(nanmean(nanmean(hcloud_mon0(:,:,8,:),4),3));
    clear hcloud_mon0
    % medium cloud
    mcloud_mon_clim_May=squeeze(nanmean(nanmean(mcloud_mon0(:,:,5,:),4),3));
    mcloud_mon_clim_Jun=squeeze(nanmean(nanmean(mcloud_mon0(:,:,6,:),4),3));
    mcloud_mon_clim_Jul=squeeze(nanmean(nanmean(mcloud_mon0(:,:,7,:),4),3));
    mcloud_mon_clim_Aug=squeeze(nanmean(nanmean(mcloud_mon0(:,:,8,:),4),3));
    clear mcloud_mon0
    % low cloud
    lcloud_mon_clim_May=squeeze(nanmean(nanmean(lcloud_mon0(:,:,5,:),4),3));
    lcloud_mon_clim_Jun=squeeze(nanmean(nanmean(lcloud_mon0(:,:,6,:),4),3));
    lcloud_mon_clim_Jul=squeeze(nanmean(nanmean(lcloud_mon0(:,:,7,:),4),3));
    lcloud_mon_clim_Aug=squeeze(nanmean(nanmean(lcloud_mon0(:,:,8,:),4),3));
    clear lcloud_mon0

    lon0(1:1360,1)=lon(81:1440,1);
    lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
    lon_ERA5=lon0; clear lon0
    lat_ERA5=lat;  clear lat

    hcloud_mon_clim_May0(   1:1360,:)=hcloud_mon_clim_May(81:1440,:);
    hcloud_mon_clim_May0(1361:1440,:)=hcloud_mon_clim_May(1:80,:); clear hcloud_mon_clim_May
    hcloud_mon_clim_May=hcloud_mon_clim_May0; clear hcloud_mon_clim_May0

    hcloud_mon_clim_Jun0(   1:1360,:)=hcloud_mon_clim_Jun(81:1440,:);
    hcloud_mon_clim_Jun0(1361:1440,:)=hcloud_mon_clim_Jun(1:80,:); clear hcloud_mon_clim_Jun
    hcloud_mon_clim_Jun=hcloud_mon_clim_Jun0; clear hcloud_mon_clim_Jun0

    hcloud_mon_clim_Jul0(   1:1360,:)=hcloud_mon_clim_Jul(81:1440,:);
    hcloud_mon_clim_Jul0(1361:1440,:)=hcloud_mon_clim_Jul(1:80,:); clear hcloud_mon_clim_Jul
    hcloud_mon_clim_Jul=hcloud_mon_clim_Jul0; clear hcloud_mon_clim_Jul0

    hcloud_mon_clim_Aug0(   1:1360,:)=hcloud_mon_clim_Aug(81:1440,:);
    hcloud_mon_clim_Aug0(1361:1440,:)=hcloud_mon_clim_Aug(1:80,:); clear hcloud_mon_clim_Aug
    hcloud_mon_clim_Aug=hcloud_mon_clim_Aug0; clear hcloud_mon_clim_Aug0


    mcloud_mon_clim_May0(   1:1360,:)=mcloud_mon_clim_May(81:1440,:);
    mcloud_mon_clim_May0(1361:1440,:)=mcloud_mon_clim_May(1:80,:); clear mcloud_mon_clim_May
    mcloud_mon_clim_May=mcloud_mon_clim_May0; clear mcloud_mon_clim_May0

    mcloud_mon_clim_Jun0(   1:1360,:)=mcloud_mon_clim_Jun(81:1440,:);
    mcloud_mon_clim_Jun0(1361:1440,:)=mcloud_mon_clim_Jun(1:80,:); clear mcloud_mon_clim_Jun
    mcloud_mon_clim_Jun=mcloud_mon_clim_Jun0; clear mcloud_mon_clim_Jun0

    mcloud_mon_clim_Jul0(   1:1360,:)=mcloud_mon_clim_Jul(81:1440,:);
    mcloud_mon_clim_Jul0(1361:1440,:)=mcloud_mon_clim_Jul(1:80,:); clear mcloud_mon_clim_Jul
    mcloud_mon_clim_Jul=mcloud_mon_clim_Jul0; clear mcloud_mon_clim_Jul0

    mcloud_mon_clim_Aug0(   1:1360,:)=mcloud_mon_clim_Aug(81:1440,:);
    mcloud_mon_clim_Aug0(1361:1440,:)=mcloud_mon_clim_Aug(1:80,:); clear mcloud_mon_clim_Aug
    mcloud_mon_clim_Aug=mcloud_mon_clim_Aug0; clear mcloud_mon_clim_Aug0


    lcloud_mon_clim_May0(   1:1360,:)=lcloud_mon_clim_May(81:1440,:);
    lcloud_mon_clim_May0(1361:1440,:)=lcloud_mon_clim_May(1:80,:); clear lcloud_mon_clim_May
    lcloud_mon_clim_May=lcloud_mon_clim_May0; clear lcloud_mon_clim_May0

    lcloud_mon_clim_Jun0(   1:1360,:)=lcloud_mon_clim_Jun(81:1440,:);
    lcloud_mon_clim_Jun0(1361:1440,:)=lcloud_mon_clim_Jun(1:80,:); clear lcloud_mon_clim_Jun
    lcloud_mon_clim_Jun=lcloud_mon_clim_Jun0; clear lcloud_mon_clim_Jun0

    lcloud_mon_clim_Jul0(   1:1360,:)=lcloud_mon_clim_Jul(81:1440,:);
    lcloud_mon_clim_Jul0(1361:1440,:)=lcloud_mon_clim_Jul(1:80,:); clear lcloud_mon_clim_Jul
    lcloud_mon_clim_Jul=lcloud_mon_clim_Jul0; clear lcloud_mon_clim_Jul0

    lcloud_mon_clim_Aug0(   1:1360,:)=lcloud_mon_clim_Aug(81:1440,:);
    lcloud_mon_clim_Aug0(1361:1440,:)=lcloud_mon_clim_Aug(1:80,:); clear lcloud_mon_clim_Aug
    lcloud_mon_clim_Aug=lcloud_mon_clim_Aug0; clear lcloud_mon_clim_Aug0
% #########################################################################


% #########################################################################
% 2023 Cloud Cover
    disp(['ERA5 Monthly HML Cloud Year# 2023'])
    load(['cloud_hml_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],'hcloud_total','mcloud_total','lcloud_total')
    hcloud_mon_2023(:,:,:)=squeeze(hcloud_total(:,1:end,:));
    mcloud_mon_2023(:,:,:)=squeeze(mcloud_total(:,1:end,:));
    lcloud_mon_2023(:,:,:)=squeeze(lcloud_total(:,1:end,:));
    clear hcloud_total mcloud_total lcloud_total
    
    
    hcloud_mon00(   1:1360,:,:)=hcloud_mon_2023(81:1440,:,:);
    hcloud_mon00(1361:1440,:,:)=hcloud_mon_2023(1:80,:,:); clear hcloud_mon_2023
    hcloud_mon_2023=hcloud_mon00; clear hcloud_mon00
    
    mcloud_mon00(   1:1360,:,:)=mcloud_mon_2023(81:1440,:,:);
    mcloud_mon00(1361:1440,:,:)=mcloud_mon_2023(1:80,:,:); clear mcloud_mon_2023
    mcloud_mon_2023=mcloud_mon00; clear mcloud_mon00
    
    lcloud_mon00(   1:1360,:,:)=lcloud_mon_2023(81:1440,:,:);
    lcloud_mon00(1361:1440,:,:)=lcloud_mon_2023(1:80,:,:); clear lcloud_mon_2023
    lcloud_mon_2023=lcloud_mon00; clear lcloud_mon00
    
    
% Cloud Cover Anomaly in 2023 from 1981-2010
    hcloud_mon_ano_May=hcloud_mon_2023(:,:,5)-hcloud_mon_clim_May;
    hcloud_mon_ano_Jun=hcloud_mon_2023(:,:,6)-hcloud_mon_clim_Jun;
    hcloud_mon_ano_Jul=hcloud_mon_2023(:,:,7)-hcloud_mon_clim_Jul;
    hcloud_mon_ano_Aug=hcloud_mon_2023(:,:,8)-hcloud_mon_clim_Aug;
    
%     hcloud_mon_2023_May=hcloud_mon_2023(:,:,5);
%     hcloud_mon_2023_Jun=hcloud_mon_2023(:,:,6);
%     hcloud_mon_2023_Jul=hcloud_mon_2023(:,:,7);
%     hcloud_mon_2023_Aug=hcloud_mon_2023(:,:,8);
    clear hcloud_mon_2023 hcloud_mon_clim*
    
    
    mcloud_mon_ano_May=mcloud_mon_2023(:,:,5)-mcloud_mon_clim_May;
    mcloud_mon_ano_Jun=mcloud_mon_2023(:,:,6)-mcloud_mon_clim_Jun;
    mcloud_mon_ano_Jul=mcloud_mon_2023(:,:,7)-mcloud_mon_clim_Jul;
    mcloud_mon_ano_Aug=mcloud_mon_2023(:,:,8)-mcloud_mon_clim_Aug;
    
%     mcloud_mon_2023_May=mcloud_mon_2023(:,:,5);
%     mcloud_mon_2023_Jun=mcloud_mon_2023(:,:,6);
%     mcloud_mon_2023_Jul=mcloud_mon_2023(:,:,7);
%     mcloud_mon_2023_Aug=mcloud_mon_2023(:,:,8);
    clear mcloud_mon_2023 mcloud_mon_clim*
    
    
    lcloud_mon_ano_May=lcloud_mon_2023(:,:,5)-lcloud_mon_clim_May;
    lcloud_mon_ano_Jun=lcloud_mon_2023(:,:,6)-lcloud_mon_clim_Jun;
    lcloud_mon_ano_Jul=lcloud_mon_2023(:,:,7)-lcloud_mon_clim_Jul;
    lcloud_mon_ano_Aug=lcloud_mon_2023(:,:,8)-lcloud_mon_clim_Aug;
    
%     lcloud_mon_2023_May=lcloud_mon_2023(:,:,5);
%     lcloud_mon_2023_Jun=lcloud_mon_2023(:,:,6);
%     lcloud_mon_2023_Jul=lcloud_mon_2023(:,:,7);
%     lcloud_mon_2023_Aug=lcloud_mon_2023(:,:,8);
    clear lcloud_mon_2023 lcloud_mon_clim*
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2. SST Climatology during 1981-2010, [\circC]
    count=1;
    for year=1981:2010
          disp(['ERA5 Monthly SST Year# ',num2str(year)])
          load(['SST_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'sst','lon','lat')
          SST0(:,:,:,count)=squeeze(sst(:,1:end,:));
          clear sst
          count=count+1;
    end
    clear year count
    SST_clim_May=squeeze(nanmean(nanmean(SST0(:,:,5,:),4),3));
    SST_clim_Jun=squeeze(nanmean(nanmean(SST0(:,:,6,:),4),3));
    SST_clim_Jul=squeeze(nanmean(nanmean(SST0(:,:,7,:),4),3));
    SST_clim_Aug=squeeze(nanmean(nanmean(SST0(:,:,8,:),4),3));
    clear SST0

    lon0(   1:1360,1)=lon(81:1440,1);
    lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
    lon_ERA5=lon0; clear lon0
    lat_ERA5=lat;  clear lat

    SST_clim_May0(   1:1360,:)=SST_clim_May(81:1440,:);
    SST_clim_May0(1361:1440,:)=SST_clim_May(1:80,:); clear SST_clim_May
    SST_clim_May=SST_clim_May0; clear SST_clim_May0

    SST_clim_Jun0(   1:1360,:)=SST_clim_Jun(81:1440,:);
    SST_clim_Jun0(1361:1440,:)=SST_clim_Jun(1:80,:); clear SST_clim_Jun
    SST_clim_Jun=SST_clim_Jun0; clear SST_clim_Jun0

    SST_clim_Jul0(   1:1360,:)=SST_clim_Jul(81:1440,:);
    SST_clim_Jul0(1361:1440,:)=SST_clim_Jul(1:80,:); clear SST_clim_Jul
    SST_clim_Jul=SST_clim_Jul0; clear SST_clim_Jul0

    SST_clim_Aug0(   1:1360,:)=SST_clim_Aug(81:1440,:);
    SST_clim_Aug0(1361:1440,:)=SST_clim_Aug(1:80,:); clear SST_clim_Aug
    SST_clim_Aug=SST_clim_Aug0; clear SST_clim_Aug0
% #########################################################################

    
% #########################################################################
% 2023 SST, [\circC]
    disp(['ERA5 Monthly SST Year# 2023'])
    load(['SST_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],'sst')
    SST0(:,:,:)=squeeze(sst(:,1:end,:));
    clear sst
    
    SST00(   1:1360,:,:)=SST0(81:1440,:,:);
    SST00(1361:1440,:,:)=SST0(1:80,:,:); clear SST0
    SST0=SST00; clear SST00
    
% SST Anomaly in 2023 from 1981-2010, [\circC]
    SST_ano_May=SST0(:,:,5)-SST_clim_May;
    SST_ano_Jun=SST0(:,:,6)-SST_clim_Jun;
    SST_ano_Jul=SST0(:,:,7)-SST_clim_Jul;
    SST_ano_Aug=SST0(:,:,8)-SST_clim_Aug;
    
    SST_2023_May=SST0(:,:,5);
    SST_2023_Jun=SST0(:,:,6);
    SST_2023_Jul=SST0(:,:,7);
    SST_2023_Aug=SST0(:,:,8);
    clear SST0 SST_clim*
% #########################################################################
% #########################################################################



%% ########################################################################
%% Mapping Low and Medium Cloud Anomaly in 2023 Relative to 1981-2010
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.165; ixe = 0.165;  ixd = 0.04; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.060; iye = 0.060;  iyd = 0.06; iyw = (1-iys-iye-1*iyd)/2;

pos{11}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('div', 'PuOr', 14,'pchip');
color0(1:12,:)=color(2:1:13,:);
color0(7,:)=(color0(7,:)+color0(8,:))./2; 
color0(5,:)=(color0(4,:)+color0(5,:))./2;    
color0(6,:)=(color0(6,:)+color0(5,:))./2;    


subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,(lcloud_mon_ano_May(:,:)'+mcloud_mon_ano_May(:,:)'));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        % #########################################################################
        % Sticks inside the SSTA contours, to indicate the warm side
        hold on
        m_plot([337 341],[48 46.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([331.2 334],[36 34.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([322 326],[25 23.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([326.5 328.5],[12 13],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([336.5 336.5],[9.2 13],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        % #########################################################################
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-0.3 0.3]);
        title('a. Low-medium cloud cover (May 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,(lcloud_mon_ano_Jun(:,:)'+mcloud_mon_ano_Jun(:,:)'));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        % #########################################################################
        % Sticks inside the SSTA contours, to indicate the warm side
        hold on
        m_plot([338 340.5],[61 58.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([322.5 325],[57.5 55.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([327 329.5],[41 44.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([329 331],[31 29],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([316 318],[22 19],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([317 317],[5 8],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([336.5 336.5],[9.2 12.2],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-0.3 0.3]);
        title('b. Low-medium cloud cover (June 2023)','fontsize',24,'FontWeight','bold')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,(lcloud_mon_ano_Jul(:,:)'+mcloud_mon_ano_Jul(:,:)'));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        % #########################################################################
        % Sticks inside the SSTA contours, to indicate the warm side
        hold on
        m_plot([330 326],[65 65],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([327.5 324.5],[55 55],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([332 329],[45 45],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([317 317],[40 43],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([307 307],[38 41],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([295 295],[39 42],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        
        hold on
        m_plot([357 360],[68 68],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([356 359],[64 64],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([352 355],[60 59],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        
        hold on
        m_plot([345 347],[45 42],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([329 331],[33 30],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([316 318],[26 23],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        
        hold on
        m_plot([317 317],[9 12],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([336.5 336.5],[10 13],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        % #########################################################################
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-0.3 0.3]);
        title('c. Low-medium cloud cover (July 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,(lcloud_mon_ano_Aug(:,:)'+mcloud_mon_ano_Aug(:,:)'));
        shading flat
        hold on
 
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        % #########################################################################
        % Sticks inside the SSTA contours, to indicate the warm side
        hold on
        m_plot([338 334],[65 65],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([331 328],[55 55],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([333 330],[45 43],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([348 345],[35 35],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([344 341],[25 25],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([320 323],[38 38],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        
        hold on
        m_plot([310 311.5],[45 48.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([365 368],[68 68],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([364 367],[65 65],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        
        hold on
        m_plot([295 295],[31 28],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([307 307],[24 21],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([316 318],[26 23],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        
        hold on
        m_plot([317 317],[9.2 12.2],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        hold on
        m_plot([336.5 336.5],[11.5 14.5],'color',[0.850, 0.325, 0.098],'linewidth',2.0,'linestyle','-')
        % #########################################################################
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-0.3 0.3]);
        title('d. Low-medium cloud cover (August 2023)','fontsize',24,'FontWeight','bold')
        
        
        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.016 iys+0.5*iyw+0*iyd 0.012 1*iyw+1*iyd]);
        set(hBar1, 'ytick',-0.3:0.1:0.3,'yticklabel',{'<-0.3','-0.2','-0.1','0','0.1','0.2','>0.3'},'fontsize',22,'FontName','Arial','LineWidth',1.2,'TickLength',0.040);

        
% #########################################################################     

disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')

