%% Extended Data Figure: Monthly MSLP in 2023 from 1981-2010 Mean
%  SLP: Sea Level Pressure, [Pa]
%  May, June, July, August, Monthly Mean 1981-2023
%  ERA-5 Monthly Averaged Reanalysis on Single Levels


%% #######################################################################
%% Figure #1: 2023 SLP Anomaly Relative to 1981-2010 Mean
clc;clear
    
% #########################################################################
% 1. MSLP Climatology during 1981-2010
    count=1;
    for year=1981:2010
          disp(['ERA5 Monthly Pressure Year# ',num2str(year)])
          load(['mslp_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'mslp','lon','lat')
          mslp0(:,:,:,count)=squeeze(mslp(:,1:end,:)./100); % Unit from [Pa] to [hPa]
          clear mslp
          count=count+1;
    end
    clear year count
    mslp_clim_May=squeeze(nanmean(nanmean(mslp0(:,:,5,:),4),3));
    mslp_clim_Jun=squeeze(nanmean(nanmean(mslp0(:,:,6,:),4),3));
    mslp_clim_Jul=squeeze(nanmean(nanmean(mslp0(:,:,7,:),4),3));
    mslp_clim_Aug=squeeze(nanmean(nanmean(mslp0(:,:,8,:),4),3));
    clear mslp0

    lon0(1:1360,1)=lon(81:1440,1);
    lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
    lon_ERA5=lon0; clear lon0
    lat_ERA5=lat;  clear lat

    mslp_clim_May0(   1:1360,:)=mslp_clim_May(81:1440,:);
    mslp_clim_May0(1361:1440,:)=mslp_clim_May(1:80,:); clear mslp_clim_May
    mslp_clim_May=mslp_clim_May0; clear mslp_clim_May0

    mslp_clim_Jun0(   1:1360,:)=mslp_clim_Jun(81:1440,:);
    mslp_clim_Jun0(1361:1440,:)=mslp_clim_Jun(1:80,:); clear mslp_clim_Jun
    mslp_clim_Jun=mslp_clim_Jun0; clear mslp_clim_Jun0

    mslp_clim_Jul0(   1:1360,:)=mslp_clim_Jul(81:1440,:);
    mslp_clim_Jul0(1361:1440,:)=mslp_clim_Jul(1:80,:); clear mslp_clim_Jul
    mslp_clim_Jul=mslp_clim_Jul0; clear mslp_clim_Jul0

    mslp_clim_Aug0(   1:1360,:)=mslp_clim_Aug(81:1440,:);
    mslp_clim_Aug0(1361:1440,:)=mslp_clim_Aug(1:80,:); clear mslp_clim_Aug
    mslp_clim_Aug=mslp_clim_Aug0; clear mslp_clim_Aug0
% #########################################################################


% #########################################################################
% 2. MSLP in 2023
    disp(['ERA5 Monthly Pressure Year# 2023'])
    load(['mslp_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_2023.mat'],'mslp')
    mslp_2023(:,:,:)=squeeze(mslp(:,1:end,:)./100); % Unit from [Pa] to [hPa]
    clear mslp

    mslp00(   1:1360,:,:)=mslp_2023(81:1440,:,:);
    mslp00(1361:1440,:,:)=mslp_2023(1:80,:,:);   clear mslp_2023
    mslp_2023=mslp00; clear mslp00
    
    
% #########################################################################
% MSLP Anomaly in 2023 from 1981-2010, [hPa]
    mslp_ano_May=mslp_2023(:,:,5)-mslp_clim_May;
    mslp_ano_Jun=mslp_2023(:,:,6)-mslp_clim_Jun;
    mslp_ano_Jul=mslp_2023(:,:,7)-mslp_clim_Jul;
    mslp_ano_Aug=mslp_2023(:,:,8)-mslp_clim_Aug;
    
    mslp_2023_May=mslp_2023(:,:,5);
    mslp_2023_Jun=mslp_2023(:,:,6);
    mslp_2023_Jul=mslp_2023(:,:,7);
    mslp_2023_Aug=mslp_2023(:,:,8);
    clear mslp_2023 mslp_clim*
% #########################################################################
% #########################################################################
    


%% ######################################################################## 
%% Mapping MSLP Anomaly in 2023 Relative to 1981-2010, [hPa]
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.165; ixe = 0.165;  ixd = 0.04; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.060; iye = 0.060;  iyd = 0.06; iyw = (1-iys-iye-1*iyd)/2;

pos{11}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:); 
color0(6,:)=(color0(6,:)+color0(6,:)+color0(5,:))./3;    


subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,mslp_ano_May');
        shading flat
        hold on
    
        % #########################################################################
        % MSLP climatology
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(mslp_2023_May,2,2)',940:2:1024,'color',[0.201, 0.645, 0.833],'linewidth',1.8);
        v21=[940:2:1060];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        m_text(345,50.5,'H','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        % #########################################################################
        
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-9 9]);
        title('a. MSLP and anomalies (May 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,mslp_ano_Jun(:,:)');
        shading flat
        hold on
        
        % #########################################################################
        % MSLP climatology
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(mslp_2023_Jun,2,2)',940:2:1060,'color',[0.201, 0.645, 0.833],'linewidth',1.8);
        v21=[940:2:1060];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        m_text(307,30,'H','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-9 9]);
        title('b. MSLP and anomalies (June 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,mslp_ano_Jul(:,:)');
        shading flat
        hold on
        
        % #########################################################################
        % MSLP climatology
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(mslp_2023_Jul,2,2)',940:2:1060,'color',[0.201, 0.645, 0.833],'linewidth',1.8);
        v21=[940:2:1060];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        m_text(328,35,'H','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        m_text(273,58,'L','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        m_text(365,61,'L','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0) 
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-9 9]);
        title('c. MSLP and anomalies (July 2023)','fontsize',24,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,mslp_ano_Aug(:,:)');
        shading flat
        hold on
        
        % #########################################################################
        % MSLP
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(mslp_2023_Aug,2,2)',940:2:1060,'color',[0.201, 0.645, 0.833],'linewidth',1.8);
        v21=[940:2:1060];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        m_text(323,38,'H','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        m_text(279,68.5,'L','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        m_text(365,61,'L','fontsize',26,'color',[.1 .1 .1],'FontWeight','bold','rotation',0)
        % #########################################################################
        
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',22,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-9 9]);
        title('d. MSLP and anomalies (August 2023)','fontsize',24,'FontWeight','bold')
        
        
        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.016 iys+0.5*iyw+0*iyd 0.012 1*iyw+1*iyd]);
        set(hBar1, 'ytick',-9.0:1.5:9.0,'yticklabel',{'<-9',[],'-6',[],'-3',[],'0',[],'3',[],'6',[],'>9'},'fontsize',22,'FontName','Arial','LineWidth',1.2,'TickLength',0.040);
        ylabel(hBar1, '[ hPa ]','rotation',90);

    
% #########################################################################     

disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')

