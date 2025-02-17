%% Main Figure 2: Surafce Wind Speed Anomaly Overlaied with Wind Arrows in 2023 Relative to 1981-2010 Mean
%  May, June, July, August, in 2023
%  ERA-5 Monthly Averaged Reanalysis on Single Levels
%  Wind Speed, U, V at 2 m


%% #######################################################################
%% Figure 1# 2023 Wind Anomaly Relative to 1981-2010 Mean
clc;clear

% #########################################################################
% 1. Wind Speed Anomaly Relative Clim Mean State
count=1;
for year=1981:2010
      disp(['loading Wind Speed... year# ',num2str(year)])
      load(['wind10m_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'uv10','lon','lat')
      uv100(:,:,:,count)=squeeze(uv10(:,1:end,:));
      clear uv10
      count=count+1;
end
clear year count
uv10_clim_May=squeeze(nanmean(nanmean(uv100(:,:,5,:),4),3));
uv10_clim_Jun=squeeze(nanmean(nanmean(uv100(:,:,6,:),4),3));
uv10_clim_Jul=squeeze(nanmean(nanmean(uv100(:,:,7,:),4),3));
uv10_clim_Aug=squeeze(nanmean(nanmean(uv100(:,:,8,:),4),3));
clear uv100

lon0(1:1360,1)=lon(81:1440,1);
lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
lon_ERA5=lon0; clear lon0
lat_ERA5=lat;  clear lat

uv10_clim_May0(1:1360,:)=uv10_clim_May(81:1440,:);
uv10_clim_May0(1361:1440,:)=uv10_clim_May(1:80,:); clear uv10_clim_May
uv10_clim_May=uv10_clim_May0; clear uv10_clim_May0

uv10_clim_Jun0(1:1360,:)=uv10_clim_Jun(81:1440,:);
uv10_clim_Jun0(1361:1440,:)=uv10_clim_Jun(1:80,:); clear uv10_clim_Jun
uv10_clim_Jun=uv10_clim_Jun0; clear uv10_clim_Jun0

uv10_clim_Jul0(1:1360,:)=uv10_clim_Jul(81:1440,:);
uv10_clim_Jul0(1361:1440,:)=uv10_clim_Jul(1:80,:); clear uv10_clim_Jul
uv10_clim_Jul=uv10_clim_Jul0; clear uv10_clim_Jul0

uv10_clim_Aug0(1:1360,:)=uv10_clim_Aug(81:1440,:);
uv10_clim_Aug0(1361:1440,:)=uv10_clim_Aug(1:80,:); clear uv10_clim_Aug
uv10_clim_Aug=uv10_clim_Aug0; clear uv10_clim_Aug0


% 2023 uv10 Anomaly relative to the clim-mean
    load(['wind10m_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],'uv10')
    uv100(:,:,:)=squeeze(uv10(:,1:end,:));
    clear uv10
    
    uv1000(1:1360,:,:)=uv100(81:1440,:,:);
    uv1000(1361:1440,:,:)=uv100(1:80,:,:); clear uv100
    uv100=uv1000; clear uv1000
    
    uv10_ano_May=uv100(:,:,5)-uv10_clim_May;
    uv10_ano_Jun=uv100(:,:,6)-uv10_clim_Jun;
    uv10_ano_Jul=uv100(:,:,7)-uv10_clim_Jul;
    uv10_ano_Aug=uv100(:,:,8)-uv10_clim_Aug;
    
    uv10_2023_May=uv100(:,:,5);
    uv10_2023_Jun=uv100(:,:,6);
    uv10_2023_Jul=uv100(:,:,7);
    uv10_2023_Aug=uv100(:,:,8);
    clear uv100 uv10_clim*
% #########################################################################

    
% #########################################################################
% 3. UV 10 m Anomaly Relative Clim Mean State
count=1;
for year=1981:2010
      disp(['loading U&V... year# ',num2str(year)])
      load(['uv10m_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'u10','v10')
      u100(:,:,:,count)=squeeze(u10(:,1:end,:));
      v100(:,:,:,count)=squeeze(v10(:,1:end,:));
      clear u10 v10
      count=count+1;
end
clear year count
u10_clim_May=squeeze(nanmean(nanmean(u100(:,:,5,:),4),3));
u10_clim_Jun=squeeze(nanmean(nanmean(u100(:,:,6,:),4),3));
u10_clim_Jul=squeeze(nanmean(nanmean(u100(:,:,7,:),4),3));
u10_clim_Aug=squeeze(nanmean(nanmean(u100(:,:,8,:),4),3));
clear u100
v10_clim_May=squeeze(nanmean(nanmean(v100(:,:,5,:),4),3));
v10_clim_Jun=squeeze(nanmean(nanmean(v100(:,:,6,:),4),3));
v10_clim_Jul=squeeze(nanmean(nanmean(v100(:,:,7,:),4),3));
v10_clim_Aug=squeeze(nanmean(nanmean(v100(:,:,8,:),4),3));
clear v100


u10_clim_May0(1:1360,:)=u10_clim_May(81:1440,:);
u10_clim_May0(1361:1440,:)=u10_clim_May(1:80,:); clear u10_clim_May
u10_clim_May=u10_clim_May0; clear u10_clim_May0

v10_clim_May0(1:1360,:)=v10_clim_May(81:1440,:);
v10_clim_May0(1361:1440,:)=v10_clim_May(1:80,:); clear v10_clim_May
v10_clim_May=v10_clim_May0; clear v10_clim_May0


u10_clim_Jun0(1:1360,:)=u10_clim_Jun(81:1440,:);
u10_clim_Jun0(1361:1440,:)=u10_clim_Jun(1:80,:); clear u10_clim_Jun
u10_clim_Jun=u10_clim_Jun0; clear u10_clim_Jun0

v10_clim_Jun0(1:1360,:)=v10_clim_Jun(81:1440,:);
v10_clim_Jun0(1361:1440,:)=v10_clim_Jun(1:80,:); clear v10_clim_Jun
v10_clim_Jun=v10_clim_Jun0; clear v10_clim_Jun0


u10_clim_Jul0(1:1360,:)=u10_clim_Jul(81:1440,:);
u10_clim_Jul0(1361:1440,:)=u10_clim_Jul(1:80,:); clear u10_clim_Jul
u10_clim_Jul=u10_clim_Jul0; clear u10_clim_Jul0

v10_clim_Jul0(1:1360,:)=v10_clim_Jul(81:1440,:);
v10_clim_Jul0(1361:1440,:)=v10_clim_Jul(1:80,:); clear v10_clim_Jul
v10_clim_Jul=v10_clim_Jul0; clear v10_clim_Jul0


u10_clim_Aug0(1:1360,:)=u10_clim_Aug(81:1440,:);
u10_clim_Aug0(1361:1440,:)=u10_clim_Aug(1:80,:); clear u10_clim_Aug
u10_clim_Aug=u10_clim_Aug0; clear u10_clim_Aug0

v10_clim_Aug0(1:1360,:)=v10_clim_Aug(81:1440,:);
v10_clim_Aug0(1361:1440,:)=v10_clim_Aug(1:80,:); clear v10_clim_Aug
v10_clim_Aug=v10_clim_Aug0; clear v10_clim_Aug0



% 2023 u10 Anomaly relative to the clim-mean
    load(['uv10m_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],'u10','v10')
    u100(:,:,:)=squeeze(u10(:,1:end,:));
    v100(:,:,:)=squeeze(v10(:,1:end,:));
    clear u10 v10
    
    u1000(1:1360,:,:)=u100(81:1440,:,:);
    u1000(1361:1440,:,:)=u100(1:80,:,:); clear u100
    u100=u1000; clear u1000
    
    v1000(1:1360,:,:)=v100(81:1440,:,:);
    v1000(1361:1440,:,:)=v100(1:80,:,:); clear v100
    v100=v1000; clear v1000
    
    u10_ano_May=u100(:,:,5)-u10_clim_May;
    u10_ano_Jun=u100(:,:,6)-u10_clim_Jun;
    u10_ano_Jul=u100(:,:,7)-u10_clim_Jul;
    u10_ano_Aug=u100(:,:,8)-u10_clim_Aug;
    
    u10_2023_May=u100(:,:,5);
    u10_2023_Jun=u100(:,:,6);
    u10_2023_Jul=u100(:,:,7);
    u10_2023_Aug=u100(:,:,8);
    clear u100 u10_clim*

    
    v10_ano_May=v100(:,:,5)-v10_clim_May;
    v10_ano_Jun=v100(:,:,6)-v10_clim_Jun;
    v10_ano_Jul=v100(:,:,7)-v10_clim_Jul;
    v10_ano_Aug=v100(:,:,8)-v10_clim_Aug;
    
    v10_2023_May=v100(:,:,5);
    v10_2023_Jun=v100(:,:,6);
    v10_2023_Jul=v100(:,:,7);
    v10_2023_Aug=v100(:,:,8);
    clear v100 v10_clim*
% #########################################################################


% #########################################################################
% 4. SST Anomaly Relative Clim Mean State
count=1;
for year=1981:2010
      disp(['loading SST... year# ',num2str(year)])
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

lon0(1:1360,1)=lon(81:1440,1);
lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
lon_ERA5=lon0; clear lon0
lat_ERA5=lat;  clear lat

SST_clim_May0(1:1360,:)=SST_clim_May(81:1440,:);
SST_clim_May0(1361:1440,:)=SST_clim_May(1:80,:); clear SST_clim_May
SST_clim_May=SST_clim_May0; clear SST_clim_May0

SST_clim_Jun0(1:1360,:)=SST_clim_Jun(81:1440,:);
SST_clim_Jun0(1361:1440,:)=SST_clim_Jun(1:80,:); clear SST_clim_Jun
SST_clim_Jun=SST_clim_Jun0; clear SST_clim_Jun0

SST_clim_Jul0(1:1360,:)=SST_clim_Jul(81:1440,:);
SST_clim_Jul0(1361:1440,:)=SST_clim_Jul(1:80,:); clear SST_clim_Jul
SST_clim_Jul=SST_clim_Jul0; clear SST_clim_Jul0

SST_clim_Aug0(1:1360,:)=SST_clim_Aug(81:1440,:);
SST_clim_Aug0(1361:1440,:)=SST_clim_Aug(1:80,:); clear SST_clim_Aug
SST_clim_Aug=SST_clim_Aug0; clear SST_clim_Aug0


% 2023 SST Anomaly relative to the clim-mean
    load(['SST_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],'sst')
    SST0(:,:,:)=squeeze(sst(:,1:end,:));
    clear sst
    
    SST00(1:1360,:,:)=SST0(81:1440,:,:);
    SST00(1361:1440,:,:)=SST0(1:80,:,:); clear SST0
    SST0=SST00; clear SST00
    
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

        


%% Plotting 2: Anomaly in 2023 Relative to 1981-2010
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.165; ixe = 0.165;  ixd = 0.04; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.060; iye = 0.060;  iyd = 0.06; iyw = (1-iys-iye-1*iyd)/2;

pos{11}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
  
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=(color(14:-1:9,:)+color(13:-1:8,:))./2;
color0(7:12,:)=color(7:-1:2,:);
color0(6,:)=(color0(6,:)+color0(5,:))./2;    


subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(uv10_ano_May(:,:)',1,1));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',4.0);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',4.0);
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
       

        % #########################################################################
        % Wind
        hold on
        u10_ano_May0=u10_ano_May;
        v10_ano_May0=v10_ano_May;
        u10_ano_May0(1325:1400,410:460)=NaN;
        v10_ano_May0(1325:1400,410:460)=NaN;
        u=smooth2a(squeeze(u10_ano_May0(1:16:end,1:16:end))*2.1,1,1);
        v=smooth2a(squeeze(v10_ano_May0(1:16:end,1:16:end))*2.1,1,1);
        lat1=lat_ERA5(1:16:end);
        lon1=lon_ERA5(1:16:end);
        [lo,la]=meshgrid(lon1,lat1); clear lon1 lat1
        coslat=cos(la.*3.14./180);
        u=u.*coslat';
        v=v.*coslat';
        
        Handle=m_quiver(lo,la,u',v',0,'k'); 
        set(Handle,'color',[.3 .3 .3],'LineWidth',1.5)
        
        hold on
        u0=2.1*5;
        v0=2.1*0;
        Handle3=m_quiver(354.5,16,u0,v0,0);
        set(Handle3, 'Color', [.1 .1 .1], 'LineWidth', 3)
        hold on
        
        m_text(353.7,20,'5 m/s','fontsize',19,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        % #########################################################################
    
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('a. Wind speed anomaly (May 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(uv10_ano_Jun(:,:)',1,1));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',4.0);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',4.0);
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
        
        
        % #########################################################################
        % Wind
        hold on
        u=smooth2a(squeeze(u10_ano_Jun(1:16:end,1:16:end))*2.1,1,1);
        v=smooth2a(squeeze(v10_ano_Jun(1:16:end,1:16:end))*2.1,1,1);
        lat1=lat_ERA5(1:16:end);
        lon1=lon_ERA5(1:16:end);
        [lo,la]=meshgrid(lon1,lat1); clear lon1 lat1
        coslat=cos(la.*3.14./180);
        u=u.*coslat';
        v=v.*coslat';
        
        Handle=m_quiver(lo,la,u',v',0,'k'); 
        set(Handle,'color',[.3 .3 .3],'LineWidth',1.5)
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. Wind speed anomaly (June 2023)','fontsize',24,'FontWeight','bold')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(uv10_ano_Jul(:,:)',1,1));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',4.0);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',4.0);
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
        
        
        % #########################################################################
        % Wind
        hold on
        u=smooth2a(squeeze(u10_ano_Jul(1:16:end,1:16:end))*2.1,1,1);
        v=smooth2a(squeeze(v10_ano_Jul(1:16:end,1:16:end))*2.1,1,1);
        lat1=lat_ERA5(1:16:end);
        lon1=lon_ERA5(1:16:end);
        [lo,la]=meshgrid(lon1,lat1); clear lon1 lat1
        coslat=cos(la.*3.14./180);
        u=u.*coslat';
        v=v.*coslat';
        
        Handle=m_quiver(lo,la,u',v',0,'k'); 
        set(Handle,'color',[.3 .3 .3],'LineWidth',1.5)
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. Wind speed anomaly (July 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(uv10_ano_Aug(:,:)',1,1));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',4.0);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',4.0);
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
       
        
        % #########################################################################
        % Wind
        hold on
        u=smooth2a(squeeze(u10_ano_Aug(1:16:end,1:16:end))*2.1,1,1);
        v=smooth2a(squeeze(v10_ano_Aug(1:16:end,1:16:end))*2.1,1,1);
        lat1=lat_ERA5(1:16:end);
        lon1=lon_ERA5(1:16:end);
        [lo,la]=meshgrid(lon1,lat1); clear lon1 lat1
        coslat=cos(la.*3.14./180);
        u=u.*coslat';
        v=v.*coslat';
        
        Handle=m_quiver(lo,la,u',v',0,'k'); 
        set(Handle,'color',[.3 .3 .3],'LineWidth',1.5)
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. Wind speed anomaly (August 2023)','fontsize',24,'FontWeight','bold')
        

        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.016 iys+0.5*iyw+0*iyd 0.012 1*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',22,'FontName','Arial','LineWidth',1.2,'TickLength',0.040);
        m_text(382,95, '[ m s^{-1} ]','fontsize',22,'FontWeight','normal')
        
        
% #########################################################################
% #########################################################################
disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')
