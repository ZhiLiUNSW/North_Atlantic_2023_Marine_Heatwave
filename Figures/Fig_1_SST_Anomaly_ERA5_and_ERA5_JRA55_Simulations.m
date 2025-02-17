%% Main Figure 1: NA SST Anomaly in 2023 Relative to 1981-2010 Mean
%  May, June, July, August, in 2023
%  ERA-5 Daily and Monthly Averaged Reanalysis on Single Levels
%  Daily SST from ACCESS-OM-2 0.25deg, ERA5 and JRA55-do runs


% #########################################################################
%% 0. Daily SST from ACCESS-OM-2 0.25deg, ERA5 and JRA55-do runs
% #########################################################################
%% a. SSTa from new ERA5 run
clc;clear

% #########################################################################
  % Basin Mask from ACCESS OM2 0.25 forced by ERA5 
    load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
        [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
        Sxy(isnan(basin_mask_NA))=NaN; 
% #########################################################################


% #########################################################################
    % Daily SST during 1981-2010
    lon_025=ncread('day_NA_sst_fields_ERA5_2023.nc','xt_ocean'); % 100w-20E
    lon_025=lon_025+360;
    lat_025=ncread('day_NA_sst_fields_ERA5_2023.nc','yt_ocean'); % 0-70N
    
    sst_daily_19812022=nan(366,length(1981:2022)');
    count=1;
    for year=1981:2022
          disp(['Daily SST Year# ',num2str(year)])
          sst_daily_NA=ncread(['day_NA_sst_fields_ERA5_',num2str(year),'.nc'],'sst'); % unit in [K]
          sst_daily_NA=sst_daily_NA-273.15; % Celsius
          
          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_19812022(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(:,1:302,:),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
          sst_daily_19812022_WNA(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(1:240,1:302,:),2),1))./squeeze(nansum(nansum(Sxy(1:240,1:302),2),1));
          sst_daily_19812022_ENA(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(241:end,1:302,:),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:302),2),1));
          clear sst_daily_NA
          count=count+1;
    end
    clear year count

    
    % 2023 Daily SST Anomaly relative to the clim-mean
          disp(['Daily SST Year# 2023'])
          sst_daily_NA=ncread(['day_NA_sst_fields_ERA5_2023.nc'],'sst'); % unit in [K]
          sst_daily_NA=sst_daily_NA-273.15; % Celsius

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_2023(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(:,1:302,1:365),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
          sst_daily_2023_WNA(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(1:240,1:302,1:365),2),1))./squeeze(nansum(nansum(Sxy(1:240,1:302),2),1));
          sst_daily_2023_ENA(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(241:end,1:302,1:365),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:302),2),1));
          sst_daily_2023(366,1)=NaN;
          sst_daily_2023_WNA(366,1)=NaN;
          sst_daily_2023_ENA(366,1)=NaN;
          clear sst_daily_NA
% #########################################################################
   

% #########################################################################
    % Anomaly is relative to the clim-daily mean
    sst_daily_19812010_clim=squeeze(nanmean(sst_daily_19812022(:,1:30),2));
    sst_daily_19812010_clim_WNA=squeeze(nanmean(sst_daily_19812022_WNA(:,1:30),2));
    sst_daily_19812010_clim_ENA=squeeze(nanmean(sst_daily_19812022_ENA(:,1:30),2));
    count=1;
    for year=1:size(sst_daily_19812022,2)
        sst_daily_19812022_ano(:,count)=sst_daily_19812022(:,count)-sst_daily_19812010_clim;
        sst_daily_19812022_ano_WNA(:,count)=sst_daily_19812022_WNA(:,count)-sst_daily_19812010_clim_WNA;
        sst_daily_19812022_ano_ENA(:,count)=sst_daily_19812022_ENA(:,count)-sst_daily_19812010_clim_ENA;
        count=count+1;
    end
    
    sst_daily_2023_ano(:,1)=sst_daily_2023(:,1)-sst_daily_19812010_clim; 
    sst_daily_2023_ano_WNA(:,1)=sst_daily_2023_WNA(:,1)-sst_daily_19812010_clim_WNA;
    sst_daily_2023_ano_ENA(:,1)=sst_daily_2023_ENA(:,1)-sst_daily_19812010_clim_ENA;
% #########################################################################


save('plot_2023_SST_Anomaly_ACCESS_OM2_025_V1_1_ERA5_SSTa_new_ERA5run.mat',...
     'sst_daily_19812022_ano','sst_daily_2023_ano','sst_daily_2023_ano_WNA','sst_daily_2023_ano_ENA')
 

 
 
% #########################################################################
%% b. SSTa from JRA5-do IAF run
 clc;clear

% #########################################################################
  % Basin Mask from ACCESS OM2 0.25 forced by ERA5
    load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
        [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
% #########################################################################


% #########################################################################
    % Daily SST during 1981-2010
    lon_025=ncread('ACCESS-OM2-025_jra55_daily_sst_field_2022.nc','xt_ocean'); % 100w-20E
    lon_025=lon_025+360;
    lat_025=ncread('ACCESS-OM2-025_jra55_daily_sst_field_2022.nc','yt_ocean'); % 0-70N
    
    sst_daily_19812022=nan(366,length(1981:2022)');
    count=1;
    for year=1981:2022
          disp(['Daily SST Year# ',num2str(year)])
          sst_daily_NA=ncread(['ACCESS-OM2-025_jra55_daily_sst_field_',num2str(year),'.nc'],'sst'); % unit in [K]
          sst_daily_NA=sst_daily_NA-273.15; % Celsius
          
          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_19812022(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(:,1:302,:),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
          sst_daily_19812022_WNA(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(1:240,1:302,:),2),1))./squeeze(nansum(nansum(Sxy(1:240,1:302),2),1));
          sst_daily_19812022_ENA(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(241:end,1:302,:),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:302),2),1));
          clear sst_daily_NA
          count=count+1;
    end
    clear year count

    
    % 2023 Daily SST Anomaly relative to the clim-mean
          disp(['Daily SST Year# 2023'])
          sst_daily_NA=ncread(['ACCESS-OM2-025_jra55_daily_sst_field_2023.nc'],'sst'); % unit in [K]
          sst_daily_NA=sst_daily_NA-273.15; % Celsius

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_2023(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(:,1:302,1:365),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
          sst_daily_2023_WNA(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(1:240,1:302,1:365),2),1))./squeeze(nansum(nansum(Sxy(1:240,1:302),2),1));
          sst_daily_2023_ENA(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(241:end,1:302,1:365),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:302),2),1));
          sst_daily_2023(366,1)=NaN;
          sst_daily_2023_WNA(366,1)=NaN;
          sst_daily_2023_ENA(366,1)=NaN;
          clear sst_daily_NA
% #########################################################################
   

% #########################################################################
    % Anomaly is relative to the clim-daily mean
    sst_daily_19812010_clim=squeeze(nanmean(sst_daily_19812022(:,1:30),2));
    sst_daily_19812010_clim_WNA=squeeze(nanmean(sst_daily_19812022_WNA(:,1:30),2));
    sst_daily_19812010_clim_ENA=squeeze(nanmean(sst_daily_19812022_ENA(:,1:30),2));
    count=1;
    for year=1:size(sst_daily_19812022,2)
        sst_daily_19812022_ano(:,count)=sst_daily_19812022(:,count)-sst_daily_19812010_clim;
        sst_daily_19812022_ano_WNA(:,count)=sst_daily_19812022_WNA(:,count)-sst_daily_19812010_clim_WNA;
        sst_daily_19812022_ano_ENA(:,count)=sst_daily_19812022_ENA(:,count)-sst_daily_19812010_clim_ENA;
        count=count+1;
    end
    
    sst_daily_2023_ano(:,1)=sst_daily_2023(:,1)-sst_daily_19812010_clim; 
    sst_daily_2023_ano_WNA(:,1)=sst_daily_2023_WNA(:,1)-sst_daily_19812010_clim_WNA;
    sst_daily_2023_ano_ENA(:,1)=sst_daily_2023_ENA(:,1)-sst_daily_19812010_clim_ENA;
% #########################################################################


save('plot_2023_SST_Anomaly_ACCESS_OM2_025_V1_2_JRA55_SSTa.mat',...
     'sst_daily_19812022_ano','sst_daily_2023_ano','sst_daily_2023_ano_WNA','sst_daily_2023_ano_ENA')
 
 


%% 1.######################################################################
%% Figure 1# Daily SST Anomaly and Map of Monthly SSTA, Relative to 1981-2010 Mean
%  Part 1: Twitter plot for daily SSTA
clc;clear

% #########################################################################
% #########################################################################
% Simulated SSTa from ACCESS OM2 0.25 forced by ERA5 and JRA55
% a. SSTa from ERA5 IAF run
load('plot_2023_SST_Anomaly_ACCESS_OM2_025_V1_1_ERA5_SSTa_new_ERA5run.mat',...
     'sst_daily_2023_ano')
sst_daily_2023_ano_ERA5=sst_daily_2023_ano; 
clear sst_daily_2023_ano
 
% b. SSTa from JRA5-do IAF run
load('plot_2023_SST_Anomaly_ACCESS_OM2_025_V1_2_JRA55_SSTa.mat',...
     'sst_daily_2023_ano')
sst_daily_2023_ano_JRA55=sst_daily_2023_ano; 
clear sst_daily_2023_ano
% #########################################################################


% #########################################################################
% SSTa from ERA5 reanalysis
% #########################################################################
  % Basin Mask from ERA5 data
    load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
    load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
       [lo,la]=meshgrid((lon)', (lat)');
       basin_mask_NA=griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       sst_daily_NA=nanmean(sst_daily_NA,3);
       basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
    
       [Sxy,~,~]=function_Cgrid_Area_Distance((lon)',(lat)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
% #########################################################################


% #########################################################################
    % Daily SST during 1981-2010
    sst_daily_19812022=nan(366,length(1981:2022)');
    count=1;
    for year=1981:2022
          disp(['Daily SST Year# ',num2str(year)])
          load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'sst_daily_NA','lon','lat')

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_19812022(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
          sst_daily_19812022_WNA(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(1:240,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(1:240,1:241),2),1));
          sst_daily_19812022_ENA(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(241:end,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:241),2),1));
          clear sst_daily_NA
          count=count+1;
    end
    clear year count

    
    % 2023 Daily SST Anomaly relative to the clim-mean
          disp(['Daily SST Year# 2023'])
          load('SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2023.mat','sst_daily')
          sst_daily_NA=sst_daily; clear sst_daily

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_2023(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,1:365),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
          sst_daily_2023_WNA(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(1:240,1:241,1:365),2),1))./squeeze(nansum(nansum(Sxy(1:240,1:241),2),1));
          sst_daily_2023_ENA(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(241:end,1:241,1:365),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:241),2),1));
          sst_daily_2023(366,1)=NaN;
          sst_daily_2023_WNA(366,1)=NaN;
          sst_daily_2023_ENA(366,1)=NaN;
          clear sst_daily_NA
% #########################################################################
   

% #########################################################################
    % Anomaly is relative to the clim-daily mean
    sst_daily_19812010_clim=squeeze(nanmean(sst_daily_19812022(:,1:30),2));
    sst_daily_19812010_clim_WNA=squeeze(nanmean(sst_daily_19812022_WNA(:,1:30),2));
    sst_daily_19812010_clim_ENA=squeeze(nanmean(sst_daily_19812022_ENA(:,1:30),2));
    count=1;
    for year=1:size(sst_daily_19812022,2)
        sst_daily_19812022_ano(:,count)=sst_daily_19812022(:,count)-sst_daily_19812010_clim;
        sst_daily_19812022_ano_WNA(:,count)=sst_daily_19812022_WNA(:,count)-sst_daily_19812010_clim_WNA;
        sst_daily_19812022_ano_ENA(:,count)=sst_daily_19812022_ENA(:,count)-sst_daily_19812010_clim_ENA;
        count=count+1;
    end
    
    sst_daily_2023_ano(:,1)=sst_daily_2023(:,1)-sst_daily_19812010_clim; 
    sst_daily_2023_ano_WNA(:,1)=sst_daily_2023_WNA(:,1)-sst_daily_19812010_clim_WNA;
    sst_daily_2023_ano_ENA(:,1)=sst_daily_2023_ENA(:,1)-sst_daily_19812010_clim_ENA;
% #########################################################################




%% ########################################################################
%  Plotting 3.1: Twitter Plot for Daily SSTA
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.300; ixe = 0.300;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.005;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

pos{101}  = [ixs          iys+1.15*iyw+1*iyd   ixw 0.85*iyw]; 
pos{102}  = [ixs          iys+0.85*iyw+1*iyd   ixw 0.286*iyw]; 

clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
subplot('position',pos{101})
    for year=1:42
        line00=plot(1:366,sst_daily_19812022_ano(:,year));
        set(line00,'color',color(year,:),'LineWidth',1.5,'linestyle','-'); 
        hold on
    end
    
    line00=plot(1:366,sst_daily_19812022_ano(:,5),'LineWidth',2);

    % JRA55-do simulated
    hold on
    line2023jra=plot(1:366,sst_daily_2023_ano_JRA55(:,1));
    set(line2023jra,'color',[0.6,0.6,0.6],'LineWidth',3.5,'linestyle','-'); 
    % ERA5 simulated
    hold on
    line2023era=plot(1:366,sst_daily_2023_ano_ERA5(:,1));
    set(line2023era,'color',[0.2,0.2,0.2],'LineWidth',3.5,'linestyle','-'); 
    hold on
    line2023=plot(1:366,sst_daily_2023_ano(:,1));
    set(line2023,'color',[0.9,0.3,0.3],'LineWidth',6,'linestyle','-'); 
    
    leg101=legend([line2023 line2023era line2023jra line00],...
               '2023','ACCESS OM2, ERA-5','ACCESS OM2, JRA-55','1981-2022',...
               'Location','northwest','NumColumns',1,'fontsize',14,'fontname','Helvetica');
    set(leg101,'fontsize',15)
    hold on
    title(leg101,'North Atlantic SST Anomaly','fontsize',16','fontname','Helvetica')
    legend('boxoff')

    
    set(gca,'Ylim',[-0.8 1.6],'ycolor','k') 
    set(gca,'YTick',-0.9:0.3:1.8)
    set(gca,'YTickLabel',{'-0.9','-0.6','-0.3','0','0.3','0.6','0.9','1.2','1.5','1.8'},'fontsize',16,'fontname','Helvetica')
    set(gca,'Xlim',[1 366]) 
    set(gca,'XTick',1:30:360)
    set(gca,'XTickLabel',[],'fontsize',16,'fontname','Helvetica')

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['Difference [ \circC ]'],'fontsize',16,'color','k','FontWeight','normal','fontname','Helvetica')

    text(-40,1.5,'a.','fontsize',18,'color','k','FontWeight','bold','fontname','Helvetica')
    
    

subplot('position',pos{102})
    line2023=plot(1:366,sst_daily_2023_ano(:,1));
    set(line2023,'color',[0.9,0.3,0.3],'LineWidth',5,'linestyle','-'); 
    hold on
    
    line_ENA=plot(1:366,sst_daily_2023_ano_ENA(:,1));
    set(line_ENA,'color',[0.9,0.3,0.3],'LineWidth',2.5,'linestyle',':'); 

    hold on
    line_WNA=plot(1:366,sst_daily_2023_ano_WNA(:,1));
    set(line_WNA,'color',[0.9,0.3,0.3],'LineWidth',2.5,'linestyle','--'); 

    leg102=legend([line2023 line_ENA line_WNA],...
               '2023','2023, Eastern NA','2023, Western NA',...
               'Location','northwest','NumColumns',1,'fontsize',14,'fontname','Helvetica');
    set(leg102,'fontsize',15,'fontname','Helvetica')
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[0.45 1.7],'ycolor','k') 
    set(gca,'YTick',0.5:0.5:1.5)
    set(gca,'YTickLabel',{'0.5','1.0','1.5'},'fontsize',16)
    set(gca,'Xlim',[1 366]) 
    set(gca,'XTick',1:30:360)
    set(gca,'XTickLabel',{'           Jan','           Feb','           Mar','           Apr','           May','           Jun',...
                          '           Jul','           Aug','           Sep','           Oct','           Nov','            Dec'},'fontsize',16,'fontname','Helvetica')
    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['[ \circC ]'],'fontsize',16,'color','k','FontWeight','normal','fontname','Helvetica')
 
    text(-40,1.6,'b.','fontsize',18,'color','k','FontWeight','bold','fontname','Helvetica')
    

    
% #########################################################################
%% Part 2: Map of 2023 Monthly SST Anomaly, May-August
% #########################################################################
% 1. SST Anomaly Relative Clim Mean State
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



%% #########################################################################  
%% Part 2: Plotting 3.2: Monthly Anomaly in 2023 Relative to 1981-2010
% #########################################################################
ixs = 0.300; ixe = 0.300;  ixd = 0.04; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.030; iye = 0.490;  iyd = 0.03; iyw = (1-iys-iye-1*iyd)/2;

pos{11}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

% #########################################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);  
color0(6,:)=(color0(6,:)+color0(5,:))./2;    
% #########################################################################


subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May(:,:)',1,1));
        shading flat
        hold on
     
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. SST anomaly (May 2023)','fontsize',16,'FontWeight','bold','fontname','Helvetica')
        
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun(:,:)',1,1));
        shading flat
        hold on
   
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. SST anomaly (June 2023)','fontsize',16,'FontWeight','bold','fontname','Helvetica')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul(:,:)',1,1));
        shading flat
        hold on
      
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. SST anomaly (July 2023)','fontsize',16,'FontWeight','bold','fontname','Helvetica')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug(:,:)',1,1));
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. SST anomaly (August 2023)','fontsize',16,'FontWeight','bold','fontname','Helvetica')
        
        
        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.014 iys+0.45*iyw+0*iyd 0.010 1.1*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',16,'fontname','Helvetica','LineWidth',1.2,'TickLength',0.058);
        m_text(385,95, '[ \circC ]','fontsize',16,'FontWeight','normal','fontname','Helvetica')
        
        
    
% #########################################################################     
disp(' Save to the built-in screen')
disp(' &')
disp(' Adjust legend text size to 14')


