%% Extended Data Figure: Monthly Air-Sea Heat Flux Anomaly in 2023 from 1981-2010 Mean
%  Qnet: Net Air-Sea Heat Fluxes; Shortwave Radiation, Latent Heat; [Watt/m2]
%  May, June, July, August, Monthly Mean 1981-2023
%  ERA-5 Monthly Averaged Reanalysis on Single Levels


%% #######################################################################
%% Figure 1# 2023 Q_net Anomaly Relative to 1981-2010 Mean
clc;clear

% #########################################################################
% #########################################################################
% 1. Net Surface Heat Flux Climatology during 1981-2010, [W/m2]
    count=1;
    for year=1981:2010
         disp(['ERA5 Monthly Air-Sea Heat Fluxes# ',num2str(year)])
         load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
               'latent','longwave','sensible','shortwave','unit_W_per_m2','lon','lat')
         heatflux(:,:,1:12,count)=latent+longwave+sensible+shortwave; % Qnet, [W/m2]
         clear latent longwave sensible shortwave
         count=count+1;
    end
    clear year count
    heatflux_clim    =squeeze(nanmean(nanmean(heatflux,4),3));
    heatflux_clim_May=squeeze(nanmean(nanmean(heatflux(:,:,5,:),4),3));
    heatflux_clim_Jun=squeeze(nanmean(nanmean(heatflux(:,:,6,:),4),3));
    heatflux_clim_Jul=squeeze(nanmean(nanmean(heatflux(:,:,7,:),4),3));
    heatflux_clim_Aug=squeeze(nanmean(nanmean(heatflux(:,:,8,:),4),3));
    clear heatflux

    lon0(   1:1360,1)=lon(81:1440,1);
    lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
    lon_ERA5=lon0; clear lon0
    lat_ERA5=lat;  clear lat

    heatflux_clim0(   1:1360,:)=heatflux_clim(81:1440,:);
    heatflux_clim0(1361:1440,:)=heatflux_clim(1:80,:); clear heatflux_clim
    heatflux_clim=heatflux_clim0; clear heatflux_clim0

    heatflux_clim_May0(   1:1360,:)=heatflux_clim_May(81:1440,:);
    heatflux_clim_May0(1361:1440,:)=heatflux_clim_May(1:80,:); clear heatflux_clim_May
    heatflux_clim_May=heatflux_clim_May0; clear heatflux_clim_May0

    heatflux_clim_Jun0(   1:1360,:)=heatflux_clim_Jun(81:1440,:);
    heatflux_clim_Jun0(1361:1440,:)=heatflux_clim_Jun(1:80,:); clear heatflux_clim_Jun
    heatflux_clim_Jun=heatflux_clim_Jun0; clear heatflux_clim_Jun0

    heatflux_clim_Jul0(   1:1360,:)=heatflux_clim_Jul(81:1440,:);
    heatflux_clim_Jul0(1361:1440,:)=heatflux_clim_Jul(1:80,:); clear heatflux_clim_Jul
    heatflux_clim_Jul=heatflux_clim_Jul0; clear heatflux_clim_Jul0

    heatflux_clim_Aug0(   1:1360,:)=heatflux_clim_Aug(81:1440,:);
    heatflux_clim_Aug0(1361:1440,:)=heatflux_clim_Aug(1:80,:); clear heatflux_clim_Aug
    heatflux_clim_Aug=heatflux_clim_Aug0; clear heatflux_clim_Aug0
% #########################################################################


% #########################################################################
% Net Surface Heat Flux in 2023, [W/m2]
    disp(['ERA5 Monthly Air-Sea Heat Fluxes# 2023'])
    load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],...
          'latent','longwave','sensible','shortwave')
    heatflux(:,:,:)=latent+longwave+sensible+shortwave; % Qnet, [W/m2]
    clear latent longwave sensible shortwave

    heatflux0(   1:1360,:,:)=heatflux(81:1440,:,:);
    heatflux0(1361:1440,:,:)=heatflux(1:80,:,:); clear heatflux
    heatflux=heatflux0; clear heatflux0
    
% #########################################################################
% Net Surface Heat Flux Anomaly in 2023 from 1981-2010, [W/m2]
    Qnet_ano_May=heatflux(:,:,5)-heatflux_clim_May;
    Qnet_ano_Jun=heatflux(:,:,6)-heatflux_clim_Jun;
    Qnet_ano_Jul=heatflux(:,:,7)-heatflux_clim_Jul;
    Qnet_ano_Aug=heatflux(:,:,8)-heatflux_clim_Aug;
    clear heatflux heatflux_clim*
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2. Shortwave Flux Climatology during 1981-2010, [W/m2]
    count=1;
    for year=1981:2010
         disp(['ERA5 Monthly Shortwave# ',num2str(year)])
         load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
               'shortwave','unit_W_per_m2','lon','lat')
         heatflux(:,:,1:12,count)=shortwave; % Qsw, [W/m2]
         clear shortwave
         count=count+1;
    end
    clear year count
    heatflux_clim    =squeeze(nanmean(nanmean(heatflux,4),3));
    heatflux_clim_May=squeeze(nanmean(nanmean(heatflux(:,:,5,:),4),3));
    heatflux_clim_Jun=squeeze(nanmean(nanmean(heatflux(:,:,6,:),4),3));
    heatflux_clim_Jul=squeeze(nanmean(nanmean(heatflux(:,:,7,:),4),3));
    heatflux_clim_Aug=squeeze(nanmean(nanmean(heatflux(:,:,8,:),4),3));
    clear heatflux

    lon0(   1:1360,1)=lon(81:1440,1);
    lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
    lon_ERA5=lon0; clear lon0
    lat_ERA5=lat;  clear lat

    heatflux_clim0(   1:1360,:)=heatflux_clim(81:1440,:);
    heatflux_clim0(1361:1440,:)=heatflux_clim(1:80,:); clear heatflux_clim
    heatflux_clim=heatflux_clim0; clear heatflux_clim0

    heatflux_clim_May0(   1:1360,:)=heatflux_clim_May(81:1440,:);
    heatflux_clim_May0(1361:1440,:)=heatflux_clim_May(1:80,:); clear heatflux_clim_May
    heatflux_clim_May=heatflux_clim_May0; clear heatflux_clim_May0

    heatflux_clim_Jun0(   1:1360,:)=heatflux_clim_Jun(81:1440,:);
    heatflux_clim_Jun0(1361:1440,:)=heatflux_clim_Jun(1:80,:); clear heatflux_clim_Jun
    heatflux_clim_Jun=heatflux_clim_Jun0; clear heatflux_clim_Jun0

    heatflux_clim_Jul0(   1:1360,:)=heatflux_clim_Jul(81:1440,:);
    heatflux_clim_Jul0(1361:1440,:)=heatflux_clim_Jul(1:80,:); clear heatflux_clim_Jul
    heatflux_clim_Jul=heatflux_clim_Jul0; clear heatflux_clim_Jul0

    heatflux_clim_Aug0(   1:1360,:)=heatflux_clim_Aug(81:1440,:);
    heatflux_clim_Aug0(1361:1440,:)=heatflux_clim_Aug(1:80,:); clear heatflux_clim_Aug
    heatflux_clim_Aug=heatflux_clim_Aug0; clear heatflux_clim_Aug0
% #########################################################################


% #########################################################################
% Shortwave Flux in 2023, [W/m2]
    disp(['ERA5 Monthly Shortwave# 2023'])
    load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],...
          'shortwave')
    heatflux(:,:,:)=shortwave; % Qsw, [W/m2]
    clear shortwave

    heatflux0(   1:1360,:,:)=heatflux(81:1440,:,:);
    heatflux0(1361:1440,:,:)=heatflux(1:80,:,:); clear heatflux
    heatflux=heatflux0; clear heatflux0
    
% #########################################################################
% Shortwave Flux Anomaly in 2023 from 1981-2010, [W/m2]
    Qsw_ano_May=heatflux(:,:,5)-heatflux_clim_May;
    Qsw_ano_Jun=heatflux(:,:,6)-heatflux_clim_Jun;
    Qsw_ano_Jul=heatflux(:,:,7)-heatflux_clim_Jul;
    Qsw_ano_Aug=heatflux(:,:,8)-heatflux_clim_Aug;
    clear heatflux heatflux_clim*
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 3. Latent Heat Flux Climatology during 1981-2010, [W/m2]
    count=1;
    for year=1981:2010
         disp(['ERA5 Monthly Latnet# ',num2str(year)])
         load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
               'latent','unit_W_per_m2','lon','lat')
         heatflux(:,:,1:12,count)=latent; % Qlat, [W/m2]
         clear latent
         count=count+1;
    end
    clear year count
    heatflux_clim    =squeeze(nanmean(nanmean(heatflux,4),3));
    heatflux_clim_May=squeeze(nanmean(nanmean(heatflux(:,:,5,:),4),3));
    heatflux_clim_Jun=squeeze(nanmean(nanmean(heatflux(:,:,6,:),4),3));
    heatflux_clim_Jul=squeeze(nanmean(nanmean(heatflux(:,:,7,:),4),3));
    heatflux_clim_Aug=squeeze(nanmean(nanmean(heatflux(:,:,8,:),4),3));
    clear heatflux

    lon0(   1:1360,1)=lon(81:1440,1);
    lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
    lon_ERA5=lon0; clear lon0
    lat_ERA5=lat;  clear lat

    heatflux_clim0(   1:1360,:)=heatflux_clim(81:1440,:);
    heatflux_clim0(1361:1440,:)=heatflux_clim(1:80,:); clear heatflux_clim
    heatflux_clim=heatflux_clim0; clear heatflux_clim0

    heatflux_clim_May0(   1:1360,:)=heatflux_clim_May(81:1440,:);
    heatflux_clim_May0(1361:1440,:)=heatflux_clim_May(1:80,:); clear heatflux_clim_May
    heatflux_clim_May=heatflux_clim_May0; clear heatflux_clim_May0

    heatflux_clim_Jun0(   1:1360,:)=heatflux_clim_Jun(81:1440,:);
    heatflux_clim_Jun0(1361:1440,:)=heatflux_clim_Jun(1:80,:); clear heatflux_clim_Jun
    heatflux_clim_Jun=heatflux_clim_Jun0; clear heatflux_clim_Jun0

    heatflux_clim_Jul0(   1:1360,:)=heatflux_clim_Jul(81:1440,:);
    heatflux_clim_Jul0(1361:1440,:)=heatflux_clim_Jul(1:80,:); clear heatflux_clim_Jul
    heatflux_clim_Jul=heatflux_clim_Jul0; clear heatflux_clim_Jul0

    heatflux_clim_Aug0(   1:1360,:)=heatflux_clim_Aug(81:1440,:);
    heatflux_clim_Aug0(1361:1440,:)=heatflux_clim_Aug(1:80,:); clear heatflux_clim_Aug
    heatflux_clim_Aug=heatflux_clim_Aug0; clear heatflux_clim_Aug0
% #########################################################################


% #########################################################################
% Latent Heat Flux in 2023, [W/m2]
    disp(['ERA5 Monthly Latent# 2023'])
    load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],...
          'latent')
    heatflux(:,:,:)=latent; % Qlat, [W/m2]
    clear latent

    heatflux0(   1:1360,:,:)=heatflux(81:1440,:,:);
    heatflux0(1361:1440,:,:)=heatflux(1:80,:,:); clear heatflux
    heatflux=heatflux0; clear heatflux0
    
% #########################################################################
% Latent Heat Flux Anomaly in 2023 from 1981-2010, [W/m2]
    Qlat_ano_May=heatflux(:,:,5)-heatflux_clim_May;
    Qlat_ano_Jun=heatflux(:,:,6)-heatflux_clim_Jun;
    Qlat_ano_Jul=heatflux(:,:,7)-heatflux_clim_Jul;
    Qlat_ano_Aug=heatflux(:,:,8)-heatflux_clim_Aug;
    clear heatflux heatflux_clim*
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 4. SST Climatology during 1981-2010, [\circC]
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

        

% #########################################################################
%% Mapping Qnet, Qsw, Qlat Anomaly in 2023 Relative to 1981-2010, [W/m2]
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.220; ixe = 0.220;  ixd = 0.035; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.035; iye = 0.035;  iyd = 0.035; iyw = (1-iys-iye-3*iyd)/4;

pos{11}  = [ixs+0*ixw+0*ixd   iys+3*iyw+3*iyd   ixw 1.0*iyw]; 
pos{12}  = [ixs+1*ixw+1*ixd   iys+3*iyw+3*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+2*ixw+2*ixd   iys+3*iyw+3*iyd   ixw 1.0*iyw]; 

pos{21}  = [ixs+0*ixw+0*ixd   iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd   iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+2*ixw+2*ixd   iys+2*iyw+2*iyd   ixw 1.0*iyw]; 

pos{31}  = [ixs+0*ixw+0*ixd   iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{32}  = [ixs+1*ixw+1*ixd   iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{33}  = [ixs+2*ixw+2*ixd   iys+1*iyw+1*iyd   ixw 1.0*iyw]; 

pos{41}  = [ixs+0*ixw+0*ixd   iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{42}  = [ixs+1*ixw+1*ixd   iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{43}  = [ixs+2*ixw+2*ixd   iys+0*iyw+0*iyd   ixw 1.0*iyw];


% #########################################################################
% Color Shading ###########################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=(color(14:-1:9,:)+color(13:-1:8,:))./2;
color0(7:12,:)=(color(6:-1:1,:)+color(7:-1:2,:))./2; 
color0(6,:)=(color0(6,:)+color0(5,:))./2;    
% #########################################################################


subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qnet_ano_May(:,:)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('a. Qnet anomaly (May 2023)','fontsize',17,'FontWeight','bold')
        
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qnet_ano_Jun(:,:)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('d. Qnet anomaly (June 2023)','fontsize',17,'FontWeight','bold')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qnet_ano_Jul(:,:)',0,0));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('g. Qnet anomaly (July 2023)','fontsize',17,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qnet_ano_Aug(:,:)',0,0));
        shading flat
        hold on
 
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('j. Qnet anomaly (August 2023)','fontsize',17,'FontWeight','bold')
% ######################################################################### 
        


% #########################################################################
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qsw_ano_May(:,:)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('b. Shortwave anomaly','fontsize',17,'FontWeight','bold')
        
        
        
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qsw_ano_Jun(:,:)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('e. Shortwave anomaly','fontsize',17,'FontWeight','bold')
        
        
subplot('position',pos{32})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qsw_ano_Jul(:,:)',0,0));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('h. Shortwave anomaly','fontsize',17,'FontWeight','bold')
        
        
        
subplot('position',pos{42})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qsw_ano_Aug(:,:)',0,0));
        shading flat
        hold on
 
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('k. Shortwave anomaly','fontsize',17,'FontWeight','bold')
% #########################################################################
        

        
% #########################################################################
subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qlat_ano_May(:,:)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('c. Latent heat anomaly','fontsize',17,'FontWeight','bold')
        
        
        
subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qlat_ano_Jun(:,:)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('f. Latent heat anomaly','fontsize',17,'FontWeight','bold')
        
        
subplot('position',pos{33})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qlat_ano_Jul(:,:)',0,0));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('i. Latent heat anomaly','fontsize',17,'FontWeight','bold')
        
        
        
subplot('position',pos{43})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(Qlat_ano_Aug(:,:)',0,0));
        shading flat
        hold on
 
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',2.5);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',2.5);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',17,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-60 60]);
        title('l. Latent heat anomaly','fontsize',17,'FontWeight','bold')
        
% #########################################################################
        
        
        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+3*ixw+2*ixd+0.014 iys+1.2*iyw+1*iyd 0.013 1.6*iyw+1*iyd]);
        set(hBar1, 'ytick',-60:10:60,'yticklabel',{'<-60',[],'-40',[],'-20',[],'0',[],'20',[],'40',[],'>60'},...
                   'fontsize',17,'FontName','Arial','LineWidth',1.2,'TickLength',0.060);
        ylabel(hBar1, '[ W m^{-2} ]','rotation',90);


% #########################################################################     

disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')

