%% Extended Data Figure: Monthly Wind Speed Anomaly in 1981-2023 Relative to 1981-2010 Monthly Climatology
%  Monthly Wind Speed 1981-2023, [m/s]
%  Anomalies are relative to the 1981-2010 monthly climatology
%  ERA-5 Monthly Averaged Reanalysis on Single Levels


%% #######################################################################
%% Figure #1: Monthly Wind Speed Anomaly from 1981 to 2023
clc;clear


% #########################################################################
% Basin Mask from Levitus et al. (2012)
    load('basin_mask_Levitus2012.mat','basin_mask','lat_Lev0','lon_Lev0')
    lon_Lev0(361,1)=380.5;
    basin_mask(361,:)=basin_mask(1,:);
    basin_mask(basin_mask>1)=NaN;
    basin_mask(basin_mask<1)=NaN;
    [lo,la]=meshgrid((260:0.25:379.75)', (0:0.25:70)');
    basin_mask_NA=griddata(lon_Lev0,lat_Lev0,basin_mask',lo',la','nearest');
    clear lo la basin_mask lat_Lev0 lon_Lev0

    [Sxy,~,~]=function_Cgrid_Area_Distance((260:0.25:379.75)',(0:0.25:70)');
    Sxy(isnan(basin_mask_NA))=NaN;
    basin_mask_NA_mon = repmat(basin_mask_NA,[1,1,12]);
% #########################################################################


% #########################################################################
% #########################################################################
% Monthly Wind Speed during 1981-2023
    wind_mon=nan(12,length(1981:2023)');
    count=1;
    for year = 1981 : 2023
          disp(['ERA5 Monthly Wind Year# ',num2str(year)])
          load(['wind10m_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'uv10','lon','lat')
          % Wind speed over the North Atlantic: 260E-20E, 0-70N (NA average is only within 0-60N)
          % North Atlantic: 260E-20E, 0-70N
          lon(1441:1520,1)=lon(1:80,1)+360; lon=lon(1041:1520,1);
          lat=lat(361:641,1);
          uv10(1441:1520,:,:)=uv10(1:80,:,:);
          uv10_mon_NA=uv10(1041:1520,361:641,:);
          clear uv10
          uv10_mon_NA(isnan(basin_mask_NA_mon))=NaN;
          
          % Area-weighted wind speed over North Atlantic: 260E-20E, 0-60N
          uv10_mon_NA=uv10_mon_NA.*Sxy;
          % NA
          wind_mon    (1:12,count)=squeeze(nansum(nansum(uv10_mon_NA(   :,   1:241,:),2),1))./squeeze(nansum(nansum(Sxy(   :,   1:241),2),1));
          % Western NA
          wind_mon_WNA(1:12,count)=squeeze(nansum(nansum(uv10_mon_NA(  1:240,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(  1:240,1:241),2),1));
          % Eastern NA
          wind_mon_ENA(1:12,count)=squeeze(nansum(nansum(uv10_mon_NA(241:end,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(241:end,1:241),2),1));
          clear uv10_mon_NA
          count=count+1;
    end
    clear year count
% #########################################################################
   

% #########################################################################
%  Monthly Wind Speed Anomaly from Monthly Climatology
    wind_mon_19812010_clim    =squeeze(nanmean(wind_mon    (:,1:30),2));
    wind_mon_19812010_clim_WNA=squeeze(nanmean(wind_mon_WNA(:,1:30),2));
    wind_mon_19812010_clim_ENA=squeeze(nanmean(wind_mon_ENA(:,1:30),2));
    count=1;
    for year=1:size(wind_mon,2)
        wind_mon_ano    (:,count)=wind_mon    (:,count)-wind_mon_19812010_clim;
        wind_mon_ano_WNA(:,count)=wind_mon_WNA(:,count)-wind_mon_19812010_clim_WNA;
        wind_mon_ano_ENA(:,count)=wind_mon_ENA(:,count)-wind_mon_19812010_clim_ENA;
        count=count+1;
    end
    
    wind_mon_2023_ano    (:,1)=wind_mon    (:,43)-wind_mon_19812010_clim;
    wind_mon_2023_ano_WNA(:,1)=wind_mon_WNA(:,43)-wind_mon_19812010_clim_WNA;
    wind_mon_2023_ano_ENA(:,1)=wind_mon_ENA(:,43)-wind_mon_19812010_clim_ENA;
% #########################################################################
% #########################################################################



%% ########################################################################
%  Plotting Monthly Wind Speed Anomaly from 1981 to 2023, in Western and Eastern NA 
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.280; ixe = 0.280;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.060; iye = 0.060;  iyd = 0.06; iyw = (1-iys-iye-1*iyd)/2;

pos{1}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{2}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
% Eastern North Atlantic Wind Speed Anomaly'
subplot('position',pos{1})
    for year=1:42
        line00=plot(1:12,wind_mon_ano_ENA(:,year));
        set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
        hold on
    end
    
    line00=plot(1:12,wind_mon_ano_ENA(:,5),'LineWidth',2);
    hold on
    line23=plot(1:12,wind_mon_2023_ano_ENA(:,1));
    set(line23,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    
    hold on
    plot(6,wind_mon_2023_ano_ENA(6,1)-0.08,'^','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2);
    hold on
    text(5.8,wind_mon_2023_ano_ENA(6,1)-0.20,'June 2023','fontsize',24,'color',[1 0.5 0],'FontWeight','normal')
    
    leg=legend([line00 line23],'1981-2022','2023','Location','northeast');
    set(leg,'fontsize',26)
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[-1 1],'ycolor','k') 
    set(gca,'YTick',-1:0.5:1)
    set(gca,'YTickLabel',{'-1.0','-0.5','0','0.5','1.0'},'fontsize',24)
    set(gca,'Xlim',[1 12]) 
    set(gca,'XTick',1:1:12)
    set(gca,'XTickLabel',[],'fontsize',24)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ m s^{-1} ]'],'fontsize',24,'color','k','FontWeight','normal')
    
    title('a. Eastern North Atlantic Wind Speed Anomaly','fontsize',24,'color','k','FontWeight','bold')


    
% Western North Atlantic Wind Speed Anomaly
subplot('position',pos{2})
    for year=1:42
        line00=plot(1:12,wind_mon_ano_WNA(:,year));
        set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
        hold on
    end
    
    line00=plot(1:12,wind_mon_ano_WNA(:,5),'LineWidth',2);
    hold on
    line23=plot(1:12,wind_mon_2023_ano_WNA(:,1));
    set(line23,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    
    hold on
    plot(7,wind_mon_2023_ano_WNA(7,1)-0.08,'^','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',12,'LineWidth',2);
    hold on
    text(6.8,wind_mon_2023_ano_WNA(7,1)-0.20,'July 2023','fontsize',24,'color',[1 0.5 0],'FontWeight','normal')
    
    leg=legend([line00 line23],'1981-2022','2023','Location','northeast');
    set(leg,'fontsize',26)
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[-1 1],'ycolor','k') 
    set(gca,'YTick',-1:0.5:1)
    set(gca,'YTickLabel',{'-1.0','-0.5','0','0.5','1.0'},'fontsize',24)
    set(gca,'Xlim',[1 12]) 
    set(gca,'XTick',1:1:12)
    set(gca,'XTickLabel',{'Jan',[],'Mar',[],'May',[],'Jul',[],'Sep',[],'Nov',[]},'fontsize',24)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ m s^{-1} ]'],'fontsize',24,'color','k','FontWeight','normal')

    title('b. Western North Atlantic Wind Speed Anomaly','fontsize',24,'color','k','FontWeight','bold')

 
% #########################################################################     

disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')


