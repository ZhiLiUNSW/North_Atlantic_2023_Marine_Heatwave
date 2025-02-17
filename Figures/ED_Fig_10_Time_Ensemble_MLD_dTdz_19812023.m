%% Extended Data Figure: Time Evolution of NA MLD and dT/dz and Maps of MLD Anomalies from Multi-Data


%% #######################################################################
%% Figure 1# 2023 MLD Anomaly Relative to 1981-2010 Mean
clc;clear
time_ann = (1980:1:2023)';
time_mon = (1980:1/12:2024-1/12)';
% #########################################################################


% #########################################################################
% NA Averaged MLD
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
    % 1. IAPv3 MLD 1980-2023
    count=1;
    for year = 1980:2023
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
         clear mld_mon Sxy_mon
         % ############################################################
         count=count+1;
    end
    clear year 
% #########################################################################



% #########################################################################
    % 2. IAPv4 MLD 1980-2023
    count=1;
    for year = 1980:2023
         disp(['Year#',num2str(year),' #Monthly MLD, IAPv4'])
         %% MLD ########################################################
         load(['mld_Rh0_0125_GSW_SA_CT_Interp_Glob_stablized_IAPv4_Temp_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP');
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
         mld_mon_NA ((count-1)*12+1:(count-1)*12+12,2) = squeeze(nansum(nansum(mld_mon(:      ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(:      ,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_NA_JJA(count,2)                       = squeeze(nansum(nansum(nanmean(mld_mon(:,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,66:125,6:8),3),2),1));
         clear mld_mon Sxy_mon
         % ############################################################
         count=count+1;
    end
    clear year 
% #########################################################################



% #########################################################################
    % 3. EN4 MLD 1980-2023
    count=1;
    for year = 1980:2023
         disp(['Year#',num2str(year),' #Monthly MLD, EN4'])
         %% MLD ########################################################
         load(['mld_Rh0_0125_GSW_SA_CT_Interp_EN4_ESM_',num2str(year),'.mat'],'mld','lon_EN4','lat_EN4')
         lat_EN4=double(lat_EN4(:,1));
         lon_EN4=double(lon_EN4(:,1));
         mld0(:,:,:)=squeeze(mld(:,:,1,1:12)); % 89.5S-89.5N
         clear mld
      
         % 64.5S - 64.5N, 20E-380E
         lon_EN40(1:340,1)=lon_EN4(21:360,1);
         lon_EN40(341:360,1)=lon_EN4(1:20,1)+360; clear lon_EN4
         lon_EN4=lon_EN40; clear lon_EN40

          mld00(1:340,:,:)=mld0(21:360,:,:);
          mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0

          % Interpolate into IAP Grids
          [lo,la]=meshgrid((21:380)', (-64.5:64.5)');
          parfor month = 1:12
              mld0=squeeze(mld00(:,:,month));
              mld_mon0=griddata(lon_EN4,lat_EN4,mld0',lo',la','linear'); 
              mld_mon(:,:,month) = mld_mon0; 
              % clear mld_mon0 mld0
          end
          clear lo la month mld00

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
         mld_mon_NA ((count-1)*12+1:(count-1)*12+12,3) = squeeze(nansum(nansum(mld_mon(:      ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(:      ,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_NA_JJA(count,3)                       = squeeze(nansum(nansum(nanmean(mld_mon(:,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,66:125,6:8),3),2),1));
         clear mld_mon Sxy_mon
         % ############################################################
         count=count+1;
    end
    clear year 
% #########################################################################



% #########################################################################
    % 4. GODAS MLD 1980-2023
    count=1;
    for year = 1980 : 2023
         disp(['Year#',num2str(year),' #Monthly MLD, GODAS'])
         %% MLD ########################################################
         mld0    = double(ncread(['dbss_obml.',num2str(year),'.nc'],'dbss_obml'));
         mld0(abs(mld0)>1e4)=NaN;
         lat_GDS = double(ncread(['dbss_obml.',num2str(year),'.nc'],'lat'));
         lon_GDS = double(ncread(['dbss_obml.',num2str(year),'.nc'],'lon'));
      
         % 64.5S - 64.5N, 20E-380E
         lon_GDS0(1:340,1)=lon_GDS(21:360,1);
         lon_GDS0(341:360,1)=lon_GDS(1:20,1)+360; clear lon_GDS
         lon_GDS=lon_GDS0; clear lon_GDS0
      
         mld00(1:340,:,:)=mld0(21:360,:,:);
         mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      
         % Interpolate into IAP Grids
         [lo,la]=meshgrid((21:380)', (-64.5:64.5)');
         parfor month = 1:12
             mld0=squeeze(mld00(:,:,month));
             mld_mon0=griddata(lon_GDS,lat_GDS,mld0',lo',la','linear'); 
             mld_mon(:,:,month) = mld_mon0; 
             % clear mld_mon0 mld0
         end
         clear lo la month mld00

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
         mld_mon_NA ((count-1)*12+1:(count-1)*12+12,4) = squeeze(nansum(nansum(mld_mon(:      ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(:      ,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_NA_JJA(count,4)                       = squeeze(nansum(nansum(nanmean(mld_mon(:,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,66:125,6:8),3),2),1));
         clear mld_mon Sxy_mon
         % ############################################################
         count=count+1;
    end
    clear year 
% #########################################################################



% #########################################################################
    % 5. ORAS5 MLD 1980-2023
    count=1;
    for year = 1980:2023
         disp(['Year#',num2str(year),' #Monthly MLD, ORAS5'])
         %% MLD ########################################################
         for month = 1:12
             if year < 2015
                 if month <= 9
                    mld0(:,:,month)  = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),'0',num2str(month),'_CONS_v0.1.nc'],'somxl030'));
                    lat_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),'0',num2str(month),'_CONS_v0.1.nc'],'nav_lat'));
                    lon_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),'0',num2str(month),'_CONS_v0.1.nc'],'nav_lon'));
                 else
                    mld0(:,:,month)  = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),num2str(month),'_CONS_v0.1.nc'],'somxl030'));
                    lat_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),num2str(month),'_CONS_v0.1.nc'],'nav_lat'));
                    lon_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),num2str(month),'_CONS_v0.1.nc'],'nav_lon'));
                 end 
             else
                 if month <= 9
                    mld0(:,:,month)  = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),'0',num2str(month),'_OPER_v0.1.nc'],'somxl030'));
                    lat_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),'0',num2str(month),'_OPER_v0.1.nc'],'nav_lat'));
                    lon_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),'0',num2str(month),'_OPER_v0.1.nc'],'nav_lon'));
                 else
                    mld0(:,:,month)  = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),num2str(month),'_OPER_v0.1.nc'],'somxl030'));
                    lat_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),num2str(month),'_OPER_v0.1.nc'],'nav_lat'));
                    lon_ORA          = double(ncread(['somxl030_control_monthly_highres_2D_',num2str(year),num2str(month),'_OPER_v0.1.nc'],'nav_lon'));
                 end 
             end
         end 
         mld0(abs(mld0)>1e4)=NaN;
         
         % 64.5S - 64.5N, 20E-380E
         lon_ORA0(   1:211,:)=lon_ORA(1232:1442,:);
         lon_ORA0( 212:639,:)=lon_ORA(   3:430, :);
         lon_ORA0(640:1440,:)=lon_ORA( 431:1231,:)+360; clear lon_ORA
         lon_ORA=lon_ORA0; clear lon_ORA0
      
         lat_ORA0(   1:211,:) =lat_ORA(1232:1442,:);
         lat_ORA0( 212:1440,:)=lat_ORA(   3:1231, :); clear lat_ORA
         lat_ORA=lat_ORA0; clear lat_ORA0


         mld00(1:211,:,:)=mld0(1232:1442,:,:);
         mld00(212:1440,:,:)=mld0(3:1231,:,:); clear mld0
      
         % Interpolate into IAP Grids
         [lo,la]=meshgrid((21:380)', (-64.5:64.5)');
         parfor month = 1:12
             mld0=squeeze(mld00(:,:,month));
             mld_mon0=griddata(lon_ORA',lat_ORA',mld0',lo',la','nearest'); 
             mld_mon(:,:,month) = mld_mon0; 
             % clear mld_mon0 mld0
         end
         clear lo la month mld00

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
         mld_mon_NA ((count-1)*12+1:(count-1)*12+12,5) = squeeze(nansum(nansum(mld_mon(:      ,66:125,:),2),1))./squeeze(nansum(nansum(Sxy_mon(:      ,66:125,:),2),1));
         % Summer Mean (JJA)
         mld_mon_NA_JJA(count,5)                       = squeeze(nansum(nansum(nanmean(mld_mon(:,66:125,6:8),3),2),1))./squeeze(nansum(nansum(nanmean(Sxy_mon(:,66:125,6:8),3),2),1));
         clear mld_mon Sxy_mon
         % ############################################################
         count=count+1;
    end
    clear year 
% #########################################################################

save R1_ED_Time_MLD_19812023_V2_IAPv3v4_EN4_GODAS_ORAS5.mat
% #########################################################################
% #########################################################################



% #########################################################################
%% Plotting ###############################################################
%  Time Evolution of JJA MLD, dT/dz from IAPv3, IAPv4, EN4, ORAS5
%  With trend lines for all products during 1980-2023
%  Map of May-Aug MLD anomay in 2023
clc; clear


% #########################################################################
% #########################################################################
% 1. Summer MLD Time Evolution from IAPv3, IAPv4, EN4, ORAS5
load R1_ED_Time_MLD_19812023_V2_IAPv3v4_EN4_GODAS_ORAS5.mat
clear lon* lat* basin* count
% #########################################################################
% #########################################################################


% #########################################################################
% #########################################################################
% 2.1 IAPv3: NA Averaged Temperature
load plot_2023_SST_Anomaly_RG18_4_Hovomuller_Temp_MLD_V1_IAPv3.mat
time_ann = (1980:1:2023)';
time_mon = (1980:1/12:2024-1/12)';


% dCT/dz [T(z=5m) ? T(z=55m) / dz(=50m)]
dTdz     = (CT_mon_NA (1,:) - CT_mon_NA (6,:))'./50; % deg-C per meter
dTdz_ENA = (CT_mon_ENA(1,:) - CT_mon_ENA(6,:))'./50; % deg-C per meter
dTdz_WNA = (CT_mon_WNA(1,:) - CT_mon_WNA(6,:))'./50; % deg-C per meter


% Annual mean and summer mean (JJA) 
count=1;
for year = 1:length(time_ann)
    dtdz_ann    (count,1) = nanmean(dTdz    ((count-1)*12+1:(count-1)*12+12),1);
    dtdz_JJA    (count,1) = nanmean(dTdz    ((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_ENA(count,1) = nanmean(dTdz_ENA((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_WNA(count,1) = nanmean(dTdz_WNA((count-1)*12+6:(count-1)*12+ 8),1);
    count=count+1;
end
clear year count dTdz* CT_mon* depth* lon* lat*
% #########################################################################


% #########################################################################
% 2.2 IAPv4: NA Averaged Temperature
load plot_2023_SST_Anomaly_RG18_4_Hovomuller_Temp_MLD_V1_IAPv4.mat
time_ann = (1980:1:2023)';
time_mon = (1980:1/12:2024-1/12)';


% dCT/dz [T(z=5m) ? T(z=55m) / dz(=50m)]
dTdz     = (CT_mon_NA (1,:) - CT_mon_NA (6,:))'./50; % deg-C per meter
dTdz_ENA = (CT_mon_ENA(1,:) - CT_mon_ENA(6,:))'./50; % deg-C per meter
dTdz_WNA = (CT_mon_WNA(1,:) - CT_mon_WNA(6,:))'./50; % deg-C per meter


% Annual mean and summer mean (JJA) 
count=1;
for year = 1:length(time_ann)
    dtdz_ann    (count,2) = nanmean(dTdz    ((count-1)*12+1:(count-1)*12+12),1);
    dtdz_JJA    (count,2) = nanmean(dTdz    ((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_ENA(count,2) = nanmean(dTdz_ENA((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_WNA(count,2) = nanmean(dTdz_WNA((count-1)*12+6:(count-1)*12+ 8),1);
    count=count+1;
end
clear year count dTdz* CT_mon* depth* lon* lat*
% #########################################################################


% #########################################################################
% 2.3 EN4-ESM: NA Averaged Temperature
load plot_2023_SST_Anomaly_RG18_4_Hovomuller_Temp_MLD_V1_EN4-ESM.mat
time_ann = (1980:1:2023)';
time_mon = (1980:1/12:2024-1/12)';


% dCT/dz [T(z=5m) ? T(z=55m) / dz(=50m)]
dTdz     = (CT_mon_NA (1,:) - CT_mon_NA (6,:))'./50; % deg-C per meter
dTdz_ENA = (CT_mon_ENA(1,:) - CT_mon_ENA(6,:))'./50; % deg-C per meter
dTdz_WNA = (CT_mon_WNA(1,:) - CT_mon_WNA(6,:))'./50; % deg-C per meter


% Annual mean and summer mean (JJA) 
count=1;
for year = 1:length(time_ann)
    dtdz_ann    (count,3) = nanmean(dTdz    ((count-1)*12+1:(count-1)*12+12),1);
    dtdz_JJA    (count,3) = nanmean(dTdz    ((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_ENA(count,3) = nanmean(dTdz_ENA((count-1)*12+6:(count-1)*12+ 8),1);
    dtdz_JJA_WNA(count,3) = nanmean(dTdz_WNA((count-1)*12+6:(count-1)*12+ 8),1);
    count=count+1;
end
clear year count dTdz* CT_mon* basin* depth* lon* lat*
% #########################################################################
% #########################################################################




% #########################################################################
% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.300; ixe = 0.300;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.005;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

pos{101}  = [ixs          iys+1.44*iyw+1*iyd   ixw 0.50*iyw]; 
pos{102}  = [ixs          iys+0.85*iyw+1*iyd   ixw 0.50*iyw]; 

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
% MLD anomaly:
mld_mon_NA_JJA = mld_mon_NA_JJA - nanmean(mld_mon_NA_JJA,1);
subplot('position',pos{101})
    line_JJA_EN4=plot(time_ann,runmean(mld_mon_NA_JJA(:,3),0));
    set(line_JJA_EN4,'color',[0.301, 0.745, 0.933],'LineWidth',3,'linestyle','-'); 
    hold on
    line_JJA_ORA=plot(time_ann,runmean(mld_mon_NA_JJA(:,5),0));
    set(line_JJA_ORA,'color',[0.929, 0.694, 0.125],'LineWidth',3,'linestyle','-'); 
    hold on
    line_JJA_IAPv4=plot(time_ann,runmean(mld_mon_NA_JJA(:,2),0));
    set(line_JJA_IAPv4,'color',[0.466, 0.674, 0.188],'LineWidth',4,'linestyle','-'); 
    hold on
    line_JJA_IAPv3=plot(time_ann,runmean(mld_mon_NA_JJA(:,1),0));
    set(line_JJA_IAPv3,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    % #########################################################################


    % #########################################################################
    % Trend Lines during 1980-2023
    time_ann=(1980:1:2023)';
    % EN4
    hold on
    mld_NA_JJA_trend = polyfit(1:length(mld_mon_NA_JJA),squeeze(mld_mon_NA_JJA(1:end,3))',1);
    y_trend_line=mld_NA_JJA_trend(2)+mld_NA_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_2023_EN4=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_2023_EN4,'color',[0.301, 0.745, 0.933],'LineWidth',2,'linestyle','--'); 
    clear y_trend_line* mld_NA_JJA_trend
    % ORAS5
    hold on
    mld_NA_JJA_trend = polyfit(1:length(mld_mon_NA_JJA),squeeze(mld_mon_NA_JJA(1:end,5))',1);
    y_trend_line=mld_NA_JJA_trend(2)+mld_NA_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_2023_ORA=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_2023_ORA,'color',[0.929, 0.694, 0.125],'LineWidth',2,'linestyle','--'); 
    clear y_trend_line* mld_NA_JJA_trend
    % IAPv4
    hold on
    mld_NA_JJA_trend = polyfit(1:length(mld_mon_NA_JJA),squeeze(mld_mon_NA_JJA(1:end,2))',1);
    y_trend_line=mld_NA_JJA_trend(2)+mld_NA_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_2023_IAPv4=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_2023_IAPv4,'color',[0.466, 0.674, 0.188],'LineWidth',2,'linestyle','--'); 
    clear y_trend_line mld_NA_JJA_trend
    % IAPv3
    hold on
    mld_NA_JJA_trend = polyfit(1:length(mld_mon_NA_JJA),squeeze(mld_mon_NA_JJA(1:end,1))',1);
    y_trend_line=mld_NA_JJA_trend(2)+mld_NA_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_2023_IAPv3=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_2023_IAPv3,'color',[0.950,0.325,0.098],'LineWidth',3,'linestyle','--'); 
    clear y_trend_line mld_NA_JJA_trend
    % #########################################################################


    % #########################################################################
    hold on
    line_JJA_IAPv3=plot(time_ann,runmean(mld_mon_NA_JJA(:,1),0));
    set(line_JJA_IAPv3,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    % #########################################################################


    leg=legend([line_JJA_IAPv3 line_JJA_IAPv4 line_JJA_ORA line_JJA_EN4 line_JJA_trend_2023_IAPv3],...
        'IAPv3','IAPv4','ORAS5 (0.03)','EN4-ESM','Linear trend over 1980-2023',...
        'Location','southwest','NumColumns',4);
    set(leg,'fontsize',14)
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[-4 4],'ycolor','k') 
    set(gca,'YTick',-4:2:5)
    set(gca,'YTickLabel',{'-4','-2','0','2','4'},'fontsize',16)
    set(gca,'Xlim',[1980-0.1 2023+0.1]) 
    set(gca,'XTick',1980:10:2020)
    set(gca,'XTickLabel',[],'fontsize',16)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['[ m ]'],'fontsize',16,'color','k','FontWeight','normal')

    title('a. Summertime (JJA) MLD anmalies','fontsize',16,'color','k','FontWeight','bold')


% Plotting ################################################################
subplot('position',pos{102})
    line_JJA_3=plot(time_ann,runmean(dtdz_JJA(:,3),0));
    set(line_JJA_3,'color',[0.301, 0.745, 0.933],'LineWidth',2,'linestyle','-'); 
    hold on
    line_JJA_2=plot(time_ann,runmean(dtdz_JJA(:,2),0));
    set(line_JJA_2,'color',[0.929, 0.694, 0.125],'LineWidth',2,'linestyle','-');  
    hold on
    line_JJA_1=plot(time_ann,runmean(dtdz_JJA(:,1),0));
    set(line_JJA_1,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    hold on
    % #########################################################################


    % #########################################################################
    % Trend Lines 1980-2023
    % IAPv3
    hold on
    dtdz_JJA_trend = polyfit(1:length(dtdz_JJA(:,1)),squeeze(dtdz_JJA(1:end,1))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=dtdz_JJA_trend(2)+dtdz_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_1=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_1,'color',[0.950,0.325,0.098],'LineWidth',3,'linestyle','--'); 
    clear y_trend_line
    % IAPv4
    hold on
    dtdz_JJA_trend = polyfit(1:length(dtdz_JJA(:,2)),squeeze(dtdz_JJA(1:end,2))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=dtdz_JJA_trend(2)+dtdz_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_2=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_2,'color',[0.929, 0.694, 0.125],'LineWidth',1.5,'linestyle','--'); 
    clear y_trend_line
    % EN4
    hold on
    dtdz_JJA_trend = polyfit(1:length(dtdz_JJA(:,3)),squeeze(dtdz_JJA(1:end,3))',1);
    time_ann=(1980:1:2023)';
    y_trend_line=dtdz_JJA_trend(2)+dtdz_JJA_trend(1)*(1:length(time_ann));  
    line_JJA_trend_3=plot(time_ann(1:end,1),y_trend_line);
    set(line_JJA_trend_3,'color',[0.301, 0.745, 0.933],'LineWidth',1.5,'linestyle','--'); 
    clear y_trend_line
    % #########################################################################


    % #########################################################################
    hold on
    line_JJA_3=plot(time_ann,runmean(dtdz_JJA(:,3),0));
    set(line_JJA_3,'color',[0.301, 0.745, 0.933],'LineWidth',3,'linestyle','-'); 
    hold on
    line_JJA_2=plot(time_ann,runmean(dtdz_JJA(:,2),0));
    set(line_JJA_2,'color',[0.929, 0.694, 0.125],'LineWidth',3,'linestyle','-');  
    hold on
    line_JJA_1=plot(time_ann,runmean(dtdz_JJA(:,1),0));
    set(line_JJA_1,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    hold on
    % #########################################################################


    leg=legend([line_JJA_1 line_JJA_2 line_JJA_3 line_JJA_trend_1],...
        'IAPv3','IAPv4','EN4-ESM','Linear trend over 1980-2023','Location','southeast','NumColumns',3);
    set(leg,'fontsize',14)
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[0.055 0.077],'ycolor','k') 
    set(gca,'YTick',0.055:0.005:0.080)
    set(gca,'YTickLabel',{'0.055','0.060','0.065','0.070','0.075','0.080'},'fontsize',16)
    set(gca,'Xlim',[1980-0.1 2023.1]) 
    set(gca,'XTick',1980:10:2020)
    set(gca,'XTickLabel',{'1980','1990','2000','2010','2020'},'fontsize',16)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ \circC m^{-1} ]'],'fontsize',16,'color','k','FontWeight','normal')

    title('b. Summertime (JJA) dT/dz anomalies, 0-50 m','fontsize',16,'color','k','FontWeight','bold')
% #########################################################################
% ######################################################################### 



% #########################################################################  
% #########################################################################
% Part 2: e-f: Changes of June-August MLD from IAPv3, IAPv4, EN4-ESM, ORAS5
load R1_ED_Map_MLD_Trends_19812023_IAPv3v4_EN4_GODAS_ORAS5.mat

mld_summer_ano = squeeze(mld_JJA_ESM(:,:,43,:)) - squeeze(nanmean(mld_JJA_ESM(:,:,1:30,:),3));


% #########################################################################
% Figure 4:
ixs = 0.300; ixe = 0.300;  ixd = 0.04; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.030; iye = 0.490;  iyd = 0.03; iyw = (1-iys-iye-1*iyd)/2;

pos{11}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=color(14:-1:9,:);
color0(7:12,:)=color(7:-1:2,:);
% #########################################################################



% #########################################################################
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 64.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_summer_ano(:,:,1)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-24 24]);
        title('c. MLD anomaly (JJA 2023, IAPv3)','fontsize',16,'FontWeight','bold')
        
        
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 64.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_summer_ano(:,:,2)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-24 24]);
        title('d. MLD anomaly (JJA 2023, IAPv4)','fontsize',16,'FontWeight','bold')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 64.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_summer_ano(:,:,3)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-24 24]);
        title('e. MLD anomaly (JJA 2023, EN4-ESM)','fontsize',16,'FontWeight','bold')
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 64.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_summer_ano(:,:,5)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-24 24]);
        title('f. MLD anomaly (JJA 2023, ORAS5)','fontsize',16,'FontWeight','bold')
        
        
        
        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.014 iys+0.45*iyw+0*iyd 0.010 1.1*iyw+1*iyd]);
        set(hBar1, 'ytick',-24:4:24,'yticklabel',{'<-24',[],'-16',[],'-8',[],'0',[],'8',[],'16',[],'>24'},'fontsize',16,'FontName','Arial','LineWidth',1.2,'TickLength',0.058);
        ylabel(hBar1, '[ m ]','rotation',90);   
    
% #########################################################################



% #########################################################################
% #########################################################################

disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')


