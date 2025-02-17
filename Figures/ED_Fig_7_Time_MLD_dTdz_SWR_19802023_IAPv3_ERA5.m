%% Extended Data Figure: Time Evolution of NA MLD, dTdz, SWR
%  Monthly MLD and T from IAP TS, 1980-2023
%  Monthly SWR from ERA5, 1980-2023
%  Summer Mean (June, July, August)

%% #######################################################################
%% Plot Time Evolution of NA MLD, dTdz, SWR, 1980-2023 - IAPv3 data
clc;clear
time_ann = (1980:1:2023)';
time_mon = (1980:1/12:2024-1/12)';
% #########################################################################


% #########################################################################
% #########################################################################
% 1. NA Averaged MLD
% #########################################################################
% Basin Mask
   load('basin_mask_Levitus2012.mat','basin_mask','lat_Lev0','lon_Lev0')
  [lo,la]=meshgrid((21:380)', (-64.5:64.5)');
   basin_mask0=griddata(lon_Lev0,lat_Lev0,basin_mask',lo',la','nearest'); clear lo la  basin_mask
   basin_mask0(basin_mask0>1)=NaN;
   basin_mask0(basin_mask0<1)=NaN;

   for month=1:12
       basin_mon(:,:,month)=basin_mask0;
   end
   clear month depth basin_mask0 lon_Lev0 lat_Lev0  
% #########################################################################
       

% #########################################################################
    % IAPv3 MLD 1980-2022
    count=1;
    for year=1980:2023
         disp(['Year#',num2str(year),' #Monthly MLD, IAPv3'])
         %% MLD ########################################################
         load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP');
         mld = squeeze(mld(:,:,1,:));
         % 64.5S - 64.5N, 20E-380E
         lat_IAP            =lat_IAP(26:155,1);

         lon0(1:340,1)      = lon_IAP(21:360,1);
         lon0(341:360,1)    = lon_IAP(1:20,1); clear lon_IAP
         lon_IAP            = lon0;            clear lon0
         lon_IAP(341:360,1) = lon_IAP(341:360,1)+360;

         mld_mon = squeeze(mld(:,26:155,1:12));  clear mld

         mld_mon0(1:340,:,:)   = mld_mon(21:360,:,:);
         mld_mon0(341:360,:,:) = mld_mon(1:20,:,:);     clear mld_mon
         mld_mon               = mld_mon0;              clear mld_mon0

         mld_mon(isnan(basin_mon))=NaN;
         % ############################################################

         % ############################################################
         %% Area
         [Sxy,~,~]      = function_Cgrid_Area_Distance(lon_IAP,lat_IAP);
         for mon=1:12
             Sxy_mon(:,:,mon)=Sxy;
         end
         clear mon Sxy
         Sxy_mon(isnan(mld_mon))   = NaN;
         Sxy_mon(isnan(basin_mon)) = NaN;

         % Volume mean MLD
         mld_mon=mld_mon.*Sxy_mon;
         % ############################################################


         % ############################################################
         %% NA Basin Averaged MLD [x y t] to [month]
         % 1.NA Averaged MLD from EQ to 60N
         mld_mon_NA ((count-1)*12+1:(count-1)*12+12,1) = squeeze(nansum(nansum(mld_mon(:      ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(:      ,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_NA_JJA(count,1)                       = squeeze(nansum(nansum(nanmean(mld_mon(:,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,66:125,6:8),3),2),1));
         % Annual Mean
         mld_mon_NA_ann(count,1)                       = squeeze(nansum(nansum(nanmean(mld_mon(:,66:125,1:12),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,66:125,1:12),3),2),1));
         % For deriving a climatology
         mld_mon_NA_clim (:,count)                     = squeeze(nansum(nansum(mld_mon(:      ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(:      ,66:125,:),2),1));
         
         % 2.Western NA
         mld_mon_WNA((count-1)*12+1:(count-1)*12+12,1) = squeeze(nansum(nansum(mld_mon(1:300  ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(1:300  ,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_WNA_JJA(count,1)                      = squeeze(nansum(nansum(nanmean(mld_mon(1:300,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(1:300,66:125,6:8),3),2),1));
         
         % 3.Eastern NA
         mld_mon_ENA((count-1)*12+1:(count-1)*12+12,1) = squeeze(nansum(nansum(mld_mon(301:end,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(301:end,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_ENA_JJA(count,1)                      = squeeze(nansum(nansum(nanmean(mld_mon(301:end,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(301:end,66:125,6:8),3),2),1));
         clear mld_mon Sxy_mon
         % ############################################################
         count=count+1;
    end
    clear year basin* count lon* lat*
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2. NA Averaged Temperature
load plot_2023_SST_Anomaly_RG18_4_Hovomuller_Temp_MLD_V1_IAPv3.mat
time_ann = (1980:1:2023)';
time_mon = (1980:1/12:2024-1/12)';


% dCT/dz [T(z=5m) ? T(z=55m) / dz(=50m)]
dTdz     = (CT_mon_NA(1,:)  - CT_mon_NA(6,:))'./50; % deg-C per meter
dTdz_ENA = (CT_mon_ENA(1,:) - CT_mon_ENA(6,:))'./50; % deg-C per meter
dTdz_WNA = (CT_mon_WNA(1,:) - CT_mon_WNA(6,:))'./50; % deg-C per meter


% Annual mean and summer mean (JJA) 
count=1;
for year = 1:length(time_ann)
    dtdz_ann(count,1) = nanmean(dTdz((count-1)*12+1:(count-1)*12+12),1);
    dtdz_JJA(count,1) = nanmean(dTdz((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_ENA(count,1) = nanmean(dTdz_ENA((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_WNA(count,1) = nanmean(dTdz_WNA((count-1)*12+6:(count-1)*12+ 8),1);
    count=count+1;
end
clear year basin* count lon* lat* CT* depth10
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 3. NA Averaged Monthly SWR from ERA5 (260-379.75E, 0-60N)
    % #########################################################################
    % Basin Mask from ERA5 data
    load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
    
    load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
       lat          = lat(1:241,1); %0-60N
       sst_daily_NA = sst_daily_NA(:,1:241,:);
       
       [lo,la]      = meshgrid((lon)', (lat)');
       basin_mask_NA= griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       sst_daily_NA = nanmean(sst_daily_NA,3);
       basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
    
       [Sxy,~,~]    = function_Cgrid_Area_Distance((lon)',(lat)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_mon      = repmat(Sxy,[1,1,12]);
    % #########################################################################


    % #########################################################################
    % Monthly SWR during 1980-2023
    count=1;
    for year = 1980 : 2023
          disp(['Monthly SWR Year# ',num2str(year),' #ERA5'])
          load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'shortwave','lon','lat')
          % 260-379.75E, 0-60N
          lon0(  1:400,1) = lon(1041:1440,1);
          lon0(401:480,1) = lon(   1:  80,1)+360;
          lon = lon0; clear lon0
          
          lat = lat(361:601,1);

          sst0(  1:400,:,:) = shortwave(1041:1440,361:601,:);
          sst0(401:480,:,:) = shortwave(   1:  80,361:601,:);
          shortwave = sst0; clear sst0
          
          shortwave(isnan(Sxy_mon))=NaN;
          shortwave_mon_NA=shortwave.*Sxy_mon; clear shortwave
          
          % NA 0-60N
          Qsw_mon_NA_JJA(count,1)=squeeze(nansum(nansum(nanmean(shortwave_mon_NA(:,1:241,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,1:241,6:8),3),2),1));
          clear shortwave_mon_NA
          count=count+1;
    end
    clear year count Sxy* lon* lat* basin* 
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
%% Plotting ###############################################################
%  Figure: JJA MLD, dt/dz, SWR in the North Atlantic
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.240; ixe = 0.240;  ixd = 0.030; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.073; iye = 0.043;  iyd = 0.045; iyw = (1-iys-iye-2*iyd)/3;

pos{11}  = [ixs     iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs     iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs     iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=color(14:-1:9,:);
color0(7:12,:)=color(7:-1:2,:);

clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);
% #########################################################################
% #########################################################################



% Plotting ################################################################
subplot('position',pos{11})
    line_JJA_WNA=plot(time_ann,runmean(mld_mon_WNA_JJA,0));
    set(line_JJA_WNA,'color',[0.301, 0.745, 0.933],'LineWidth',2,'linestyle','-'); 
    hold on
    line_JJA_ENA=plot(time_ann,runmean(mld_mon_ENA_JJA,0));
    set(line_JJA_ENA,'color',[0.929, 0.694, 0.125],'LineWidth',2,'linestyle','-'); 
    hold on
    line_JJA=plot(time_ann,runmean(mld_mon_NA_JJA,0));
    set(line_JJA,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    % #########################################################################

    % NA MLD within 1980-2000, 1990-2000
    mld_mon_NA_JJA_19801990 = nanmean(mld_mon_NA_JJA( 1:11,1),1);
    mld_mon_NA_JJA_19802000 = nanmean(mld_mon_NA_JJA( 1:21,1),1);
    mld_mon_NA_JJA_19902000 = nanmean(mld_mon_NA_JJA(11:21,1),1);
    mld_mon_NA_JJA_20002010 = nanmean(mld_mon_NA_JJA(21:31,1),1);
    mld_mon_NA_JJA_20102020 = nanmean(mld_mon_NA_JJA(31:41,1),1);
    %
    
    % #########################################################################
    % Trend Lines
    % 2020
    hold on
    mld_NA_JJA_trend = polyfit(1:length(mld_mon_NA_JJA)-3,squeeze(mld_mon_NA_JJA(1:end-3,1))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=mld_NA_JJA_trend(2)+mld_NA_JJA_trend(1)*(1:length(time_ann)-3);  
    line_JJA_trend_2020=plot(time_ann(1:end-3,1),y_trend_line);
    set(line_JJA_trend_2020,'color',[0.8,0.8,0.8],'LineWidth',6,'linestyle','-'); 
    clear y_trend_line
    
    % 2022
    hold on
    mld_NA_JJA_trend = polyfit(1:length(mld_mon_NA_JJA)-1,squeeze(mld_mon_NA_JJA(1:end-1,1))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=mld_NA_JJA_trend(2)+mld_NA_JJA_trend(1)*(1:length(time_ann)-1);  
    line_JJA_trend_2022=plot(time_ann(1:end-1,1),y_trend_line);
    set(line_JJA_trend_2022,'color',[0.950,0.325,0.098],'LineWidth',3,'linestyle','--'); 
    clear y_trend_line
    % #########################################################################


    % #########################################################################
    hold on
    line_JJA=plot(time_ann,runmean(mld_mon_NA_JJA,0));
    set(line_JJA,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    % #########################################################################


    leg=legend([line_JJA line_JJA_WNA line_JJA_ENA line_JJA_trend_2020 line_JJA_trend_2022],'Summer NA MLD',...
        'Summer western NA MLD','Summer eastern NA MLD','Linear trend over 1980-2020','Linear trend over 1980-2022','Location','southwest','NumColumns',3);
    set(leg,'fontsize',16)
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[17.5 30.5],'ycolor','k') 
    set(gca,'YTick',15:5:30)
    set(gca,'YTickLabel',{'15','20','25','30'},'fontsize',19)
    set(gca,'Xlim',[1980-0.1 2023+0.1]) 
    set(gca,'XTick',1980:10:2020)
    set(gca,'XTickLabel',[],'fontsize',19)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['[ m ]'],'fontsize',20,'color','k','FontWeight','normal')

    title('a. Summertime (JJA) MLD, 1980-2023, North Atlantic','fontsize',20,'color','k','FontWeight','bold')

        

% Plotting ################################################################
subplot('position',pos{21})
    line_JJA_WNA=plot(time_ann,runmean(dtdz_JJA_WNA,0));
    set(line_JJA_WNA,'color',[0.301, 0.745, 0.933],'LineWidth',2,'linestyle','-'); 
    hold on
    line_JJA_ENA=plot(time_ann,runmean(dtdz_JJA_ENA,0));
    set(line_JJA_ENA,'color',[0.929, 0.694, 0.125],'LineWidth',2,'linestyle','-');  
    hold on
    line_JJA=plot(time_ann,runmean(dtdz_JJA,0));
    set(line_JJA,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    hold on
    % #########################################################################


    % #########################################################################
    % Trend Lines
    % 2020
    hold on
    dtdz_JJA_trend = polyfit(1:length(dtdz_JJA)-3,squeeze(dtdz_JJA(1:end-3,1))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=dtdz_JJA_trend(2)+dtdz_JJA_trend(1)*(1:length(time_ann)-3);  
    line_JJA_trend_2020=plot(time_ann(1:end-3,1),y_trend_line);
    set(line_JJA_trend_2020,'color',[0.8,0.8,0.8],'LineWidth',6,'linestyle','-'); 
    clear y_trend_line
        
    % 2022
    hold on
    dtdz_JJA_trend = polyfit(1:length(dtdz_JJA)-1,squeeze(dtdz_JJA(1:end-1,1))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=dtdz_JJA_trend(2)+dtdz_JJA_trend(1)*(1:length(time_ann)-1);  
    line_JJA_trend_2022=plot(time_ann(1:end-1,1),y_trend_line);
    set(line_JJA_trend_2022,'color',[0.950,0.325,0.098],'LineWidth',3,'linestyle','--'); 
    clear y_trend_line
    
    text(1994,0.058,'0.002 \circC m^{-1} per decade','fontsize',16,'color',[0.950,0.325,0.098],'FontWeight','normal')
    % #########################################################################


    % #########################################################################
    hold on
    line_JJA=plot(time_ann,runmean(dtdz_JJA,0));
    set(line_JJA,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    hold on
    % #########################################################################


    leg=legend([line_JJA line_JJA_WNA line_JJA_ENA line_JJA_trend_2020 line_JJA_trend_2022],'Summer NA dT/dz',...
        'Summer western NA dT/dz','Summer eastern NA dT/dz','Linear trend over 1980-2020','Linear trend over 1980-2022','Location','northwest','NumColumns',3);
    set(leg,'fontsize',16)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')

    
    set(gca,'Ylim',[0.055 0.080],'ycolor','k') 
    set(gca,'YTick',0.055:0.005:0.080)
    set(gca,'YTickLabel',{'0.055','0.060','0.065','0.070','0.075','0.080'},'fontsize',19)
    set(gca,'Xlim',[1980-0.1 2023.1]) 
    set(gca,'XTick',1980:10:2020)
    set(gca,'XTickLabel',[],'fontsize',19)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['[ \circC m^{-1} ]'],'fontsize',20,'color','k','FontWeight','normal')

    title('b. Summertime (JJA) dT/dz, 1980-2023, North Atlantic','fontsize',20,'color','k','FontWeight','bold')
% #########################################################################

    

% #########################################################################
% Plotting ################################################################
% SWR anomaly:
Qsw_mon_NA_JJA = Qsw_mon_NA_JJA - nanmean(Qsw_mon_NA_JJA(1:end,:),1);

subplot('position',pos{31})
    line_0=plot(time_ann,ones(length(time_ann))*0);
    set(line_0,'color',[0.5, 0.5, 0.5],'LineWidth',1.5,'linestyle','-'); 
    hold on

    line_JJA_1=plot(time_ann,runmean(Qsw_mon_NA_JJA(:,1),0));
    set(line_JJA_1,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    % #########################################################################


    % #########################################################################
    % Trend Lines for 1980-2023
    % ERA5
    hold on
    Qsw_NA_JJA_trend = polyfit(1:size(Qsw_mon_NA_JJA,1),squeeze(Qsw_mon_NA_JJA(1:end,1))',1);
    disp(['Qsw trend#',num2str(Qsw_NA_JJA_trend(1)*10),' W/m2 per decade'])
    time_ann=(1980:1:2023)';
    y_trend_line=Qsw_NA_JJA_trend(2)+Qsw_NA_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_1=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_1,'color',[0.950,0.325,0.098],'LineWidth',3,'linestyle','--'); 
    clear y_trend_line
    
    text(1994,-2,'0.58 W m^{-2} per decade','fontsize',16,'color',[0.950,0.325,0.098],'FontWeight','normal')
    % #########################################################################

    
    % #########################################################################
    hold on
    line_JJA_1=plot(time_ann,runmean(Qsw_mon_NA_JJA(:,1),0));
    set(line_JJA_1,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    % #########################################################################

    
    leg=legend([line_JJA_1 line_JJA_trend_1 ],...
        'Summer NA Qsw','Linear trend over 1980-2023',...
        'Location','northwest','NumColumns',2);
    set(leg,'fontsize',16)
    hold on
    % title(leg,'Jun-Aug Qsw anomaly','fontsize',17')
    legend('boxoff')
    
    set(gca,'Ylim',[-6 6],'ycolor','k') 
    set(gca,'YTick',-6:2:6)
    set(gca,'YTickLabel',{'-6','-4','-2','0','2','4','6'},'fontsize',19)
    set(gca,'Xlim',[1980-0.1 2023+0.1]) 
    set(gca,'XTick',1980:10:2020)
    set(gca,'XTickLabel',{'1980','1990','2000','2010','2020'},'fontsize',19)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['[ W m^{-2} ]'],'fontsize',20,'color','k','FontWeight','normal')
    xlabel(['Year'],'fontsize',20,'color','k','FontWeight','normal')
    
    title('c. Summertime (JJA) surface shortwave, 1980-2023, North Atlantic','fontsize',20,'color','k','FontWeight','bold')

    
        
% #########################################################################   
% #########################################################################

disp(' Save to the built-in screen')
disp(' >>')




