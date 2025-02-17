%% Main Figure 5: MLT Decomposition from IAP and ERA5
%  1981-2023:
%  Monthly MLT from IAP and Air-Sea Fluxes from ERA5
%  Interpolate monthly heat flux data to 1 degree grids
%  Monthly heat flux climatology during 1981-2010


%% ########################################################################
%  ########################################################################
% %% 1. Monthly MLT from IAP for 1981-2023
% %  Monthly MLT 1981-2023 
% clc;clear
% time_ann=(1981:2023)';
% 
% % #########################################################################
% count_yr=1;
% for year=1981:2023
%     % #####################################################################
%     % Mixed-Layer Temp from IAP
%     % Monthly MLD from IAP
%     disp(['   MLT from IAPv3, Year#',num2str(year)])
%       % dRh0=0.125 ########################################################
%       load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%       lat_IAP=lat_IAP(:,1);
%       mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%       clear mld
%       
%       lon_IAP0(1:340,1)=lon_IAP(21:360,1);
%       lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
%       lon_IAP=lon_IAP0; clear lon_IAP0
%       
%       mld00(1:340,:,:)=mld0(21:360,:,:);
%       mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
%       mld_monthly=mld00; clear mld00
%     % #####################################################################
%       
%       
%     % #####################################################################
%     % Monthly MLT from IAP
%       load(['CT_depth_Monthly_gswTSinterp_Glob_10itvl_IAP_C17_',num2str(year),'.mat'],'CT_depth','depth10')
%       CT_depth0(1:340,1:180,1:196,1:12)=CT_depth(21:360,1:180,1:196,1:12);
%       CT_depth0(341:360,1:180,1:196,1:12)=CT_depth(1:20,1:180,1:196,1:12);
%       clear CT_depth
%       depth10=depth10(1:196,1);
% 
%       depth_MML=nan(length(lon_IAP),length(lat_IAP),length(depth10),12);
%       for month=1:12
%           for i=1:length(lon_IAP)
%               for j=1:length(lat_IAP)
%                   k_MML=find(abs(depth10(:,1)-mld_monthly(i,j,month))==min(abs(depth10(:,1)-mld_monthly(i,j,month))));
%                   depth_MML(i,j,1:k_MML,month)=1;
%                   clear k_MML
%               end
%           end
%       end
%       clear i j month
% 
%       CT_depth0(isnan(depth_MML))=NaN;
%       mlt_monthly(:,:,:)=squeeze(nanmean(CT_depth0,3));
%       clear CT_depth0 
%       clear depth_MML depth10 mld_monthly mld_10m
%     % #####################################################################
% 
%       save(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],...
%             'mlt_monthly','lon_IAP','lat_IAP')
%       clear  mlt_monthly
%       count_yr=count_yr+1;
% end
% % #########################################################################
% % #########################################################################
% 
% 
% 
% %% #########################################################################
% %  #########################################################################
% % 2.1 MLT Budget in 1981-2023
% clc;clear
% time_ann=(1981:2010)';
% 
% % #########################################################################
%   % ACCESS Om2 0.25
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon lat
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% % #########################################################################
%       % second_2_month=60*60*24*30;
%       second_2_month(1,1)  = 60*60*24*31;
%       second_2_month(2,1)  = 60*60*24*28.25;
%       second_2_month(3,1)  = 60*60*24*31;
%       second_2_month(4,1)  = 60*60*24*30;
%       second_2_month(5,1)  = 60*60*24*31;
%       second_2_month(6,1)  = 60*60*24*30;
%       second_2_month(7,1)  = 60*60*24*31;
%       second_2_month(8,1)  = 60*60*24*31;
%       second_2_month(9,1)  = 60*60*24*30;
%       second_2_month(10,1) = 60*60*24*31;
%       second_2_month(11,1) = 60*60*24*30;
%       second_2_month(12,1) = 60*60*24*31;
%       
%       for month=1:12
%           second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
%       end
%       second_2_month=second_2_month0; clear second_2_month0
% % #########################################################################
%       
% 
% % #########################################################################
% count_yr=1;
% for year=1981:2023
%     disp([' Predicting NA daily SST Year#',num2str(year)])
%     % #####################################################################
%     % Monthly MLD from IAP
%     disp(['   MLD from IAPv3, Year#',num2str(year)])
%       % dRh0=0.125 ########################################################
%       load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%       lat_IAP=lat_IAP(:,1);
%       mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%       clear mld
%       
%       lon_IAP0(1:340,1)  =lon_IAP(21:360,1);
%       lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
%       lon_IAP            =lon_IAP0;            clear lon_IAP0
%       
%       mld00(1:340,:,:)   =mld0(21:360,:,:);
%       mld00(341:360,:,:) =mld0(1:20,:,:); clear mld0
%       mld_monthly        =mld00;          clear mld00
%     % #####################################################################
%         
%     % Monthly heat flux terms
%       load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
%             'shortwave','latent','longwave','sensible') 
%       Qnet=latent+longwave+sensible+shortwave;
% 
%       % ##################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_monthly;
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld = - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
%         Qnet      = Qnet      + shortwave_mld;
%         shortwave = shortwave + shortwave_mld;
%         clear shortwave_mld
%       % ##################################################################   
%        
%        
%       % ##################################################################
%       % Mixed layer warming due to Qnet
%         disp(['   Projected ML Warming in Year#',num2str(year)])
%         % Heat flux term in heat budget
%         MLT_Qnet_mon(:,:,:)  = Qnet     ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear Qnet
%         MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear shortwave
%         MLT_Qlat_mon(:,:,:)  = latent   ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear latent
%         MLT_Qlon_mon(:,:,:)  = longwave ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear longwave
%         MLT_Qsen_mon(:,:,:)  = sensible ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear sensible
%         clear mld_monthly
% 
%         
%          % ################################################################
%          % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
%          dMLT_Qnet_mon0=MLT_Qnet_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qswr_mon0=MLT_Qswr_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlat_mon0=MLT_Qlat_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlon_mon0=MLT_Qlon_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qsen_mon0=MLT_Qsen_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          Sxy_NA0=Sxy_NA;
%          for month=2:12
%             Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
%          end
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          dMLT_Qnet_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlat_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlon_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qsen_mon0(isnan(Sxy_NA0))=NaN; 
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          
%          dMLT_Qnet_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qswr_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlat_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlat_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlon_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlon_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qsen_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qsen_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          clear dMLT_*_mon0 Sxy_NA0
%       % ###################################################################
%      
%          save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'.mat'],...
%               'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
%               'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%          clear dMLT_*_mon
%          % ################################################################
%       count_yr=count_yr+1;
% end
% clear *_dail_clim
% % #########################################################################
% % #########################################################################  
% 
% 
% 
% % #########################################################################
% % #########################################################################  
% %% 2.2 MLT Budget in 1981-2023 - MLD' Term, Using SWR Climatology Everywhere
% clc;clear
% time_ann=(1981:2010)';
% 
% % #########################################################################
%   % ACCESS Om2 0.25
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon lat
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% 
% % #########################################################################
%       % second_2_month=60*60*24*30;
%       second_2_month(1,1)  = 60*60*24*31;
%       second_2_month(2,1)  = 60*60*24*28.25;
%       second_2_month(3,1)  = 60*60*24*31;
%       second_2_month(4,1)  = 60*60*24*30;
%       second_2_month(5,1)  = 60*60*24*31;
%       second_2_month(6,1)  = 60*60*24*30;
%       second_2_month(7,1)  = 60*60*24*31;
%       second_2_month(8,1)  = 60*60*24*31;
%       second_2_month(9,1)  = 60*60*24*30;
%       second_2_month(10,1) = 60*60*24*31;
%       second_2_month(11,1) = 60*60*24*30;
%       second_2_month(12,1) = 60*60*24*31;
%       
%       for month=1:12
%           second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
%       end
%       second_2_month=second_2_month0; clear second_2_month0
% % #########################################################################
% 
% 
% % #########################################################################
% count_yr=1;
% for year=1981:2023
%     disp([' Predicting NA daily SST Year#',num2str(year)])
%     % #####################################################################
%     % Monthly MLD from IAP
%     disp(['   MLD from IAPv3, Year#',num2str(year)])
%       % dRh0=0.125 ########################################################
%       load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%       lat_IAP=lat_IAP(:,1);
%       mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%       clear mld
%       
%       lon_IAP0(1:340,1)   = lon_IAP(21:360,1);
%       lon_IAP0(341:360,1) = lon_IAP(1:20,1)+360; clear lon_IAP
%       lon_IAP             = lon_IAP0;            clear lon_IAP0
%       
%       mld00(1:340,:,:)    = mld0(21:360,:,:);
%       mld00(341:360,:,:)  = mld0(1:20,:,:); clear mld0
%       mld_monthly         = mld00;          clear mld00
%     % #####################################################################
%     
%        
%     % #####################################################################
%     % Qnet and SWR Climatology
%       count=1;
%       for year_qnet=1981:2010
%           disp(['    heat fluxes Year#',num2str(year_qnet)])
%           load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year_qnet),'_1deg.mat'],...
%                 'shortwave','latent','longwave','sensible') 
%           shortwave0(:,:,:,count) = shortwave;
%           latent0   (:,:,:,count) = latent;
%           longwave0 (:,:,:,count) = longwave;
%           sensible0 (:,:,:,count) = sensible;
%           Qnet0     (:,:,:,count) = latent+longwave+sensible+shortwave;
%           clear latent longwave sensible shortwave
%           count=count+1;
%       end
%       clear year_qnet count
%       Qnet      = nanmean(Qnet0,     4); clear Qnet0
%       shortwave = nanmean(shortwave0,4); clear shortwave0
%       latent    = nanmean(latent0,   4); clear latent0
%       longwave  = nanmean(longwave0, 4); clear longwave0
%       sensible  = nanmean(sensible0, 4); clear sensible0
% 
%       % ##################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_monthly;
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld = - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
%         Qnet      = Qnet      + shortwave_mld;
%         shortwave = shortwave + shortwave_mld;
%         clear shortwave_mld
%       % ##################################################################   
%        
%        
%       % ##################################################################
%       % Mixed layer warming due to Qnet
%         disp(['   Projected ML Warming in Year#',num2str(year)])
%         % Heat flux term in heat budget
%         MLT_Qnet_mon(:,:,:)  = Qnet     ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear Qnet
%         MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear shortwave
%         MLT_Qlat_mon(:,:,:)  = latent   ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear latent
%         MLT_Qlon_mon(:,:,:)  = longwave ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear longwave
%         MLT_Qsen_mon(:,:,:)  = sensible ./mld_monthly./3992./1027.*second_2_month;  % K/month
%         clear sensible
%         clear mld_monthly
%         
%          % ################################################################
%          % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
%          dMLT_Qnet_mon0=MLT_Qnet_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qswr_mon0=MLT_Qswr_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlat_mon0=MLT_Qlat_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlon_mon0=MLT_Qlon_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qsen_mon0=MLT_Qsen_mon(240:360,91:150,:).*Sxy_NA(:,:);
% 
%          Sxy_NA0=Sxy_NA;
%          for month=2:12
%             Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
%          end
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          dMLT_Qnet_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlat_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlon_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qsen_mon0(isnan(Sxy_NA0))=NaN; 
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          % figure;imagesc(dMLT_Qnet_mon0(:,:,1)); 
%          % figure;imagesc(Sxy_NA0(:,:,1)); 
%          
%          dMLT_Qnet_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qswr_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlat_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlat_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlon_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlon_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qsen_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qsen_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          clear dMLT_*_mon0 Sxy_NA0
%       % ###################################################################
% 
%          save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'_2_MLD_prime.mat'],...
%               'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
%               'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%           clear NA_SST_mon* dMLT_*_mon
%          % ################################################################
%       count_yr=count_yr+1;
% end
% clear *_dail_clim
% % #########################################################################
% % #########################################################################  
% 
% 
% 
% 
% %% 2.3 MLT Budget in 1981-2023 - Qsw' Term, Using MLD Climatology Everywhere
% clc;clear
% time_ann=(1981:2010)';
% 
% % #########################################################################
%   % ACCESS OM2 0.25
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon lat
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% 
% % #########################################################################
%       % second_2_month=60*60*24*30;
%       second_2_month(1,1)  = 60*60*24*31;
%       second_2_month(2,1)  = 60*60*24*28.25;
%       second_2_month(3,1)  = 60*60*24*31;
%       second_2_month(4,1)  = 60*60*24*30;
%       second_2_month(5,1)  = 60*60*24*31;
%       second_2_month(6,1)  = 60*60*24*30;
%       second_2_month(7,1)  = 60*60*24*31;
%       second_2_month(8,1)  = 60*60*24*31;
%       second_2_month(9,1)  = 60*60*24*30;
%       second_2_month(10,1) = 60*60*24*31;
%       second_2_month(11,1) = 60*60*24*30;
%       second_2_month(12,1) = 60*60*24*31;
%       
%       for month=1:12
%           second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
%       end
%       second_2_month=second_2_month0; clear second_2_month0 month
% % #########################################################################
%       
% 
% % #########################################################################
%   % Monthly MLD Climatology from IAP
%     count=1;
%     for year=1981:2010
%         disp(['   MLD from IAPv3, Year#',num2str(year)])
%           % dRh0=0.125 ########################################################
%           load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%           lat_IAP=lat_IAP(:,1);
%           mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%           clear mld
% 
%           lon_IAP0(1:340,1)   = lon_IAP(21:360,1);
%           lon_IAP0(341:360,1) = lon_IAP(1:20,1)+360; clear lon_IAP
%           lon_IAP             = lon_IAP0;            clear lon_IAP0
% 
%           mld00(1:340,:,:)    = mld0(21:360,:,:);
%           mld00(341:360,:,:)  = mld0(1:20,:,:); clear mld0
%           mld_monthly(:,:,:,count)=mld00;       clear mld00
%           count=count+1;
%     end
%     clear year count
%     mld_clim=nanmean(mld_monthly,4); clear mld_monthly
% % #########################################################################
% 
% 
% % #########################################################################
% count_yr=1;
% for year=1981:2023
%     disp([' Predicting NA daily SST Year#',num2str(year)])
%     % ##################################################################### 
%     % Monthly heat flux terms
%       load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
%             'shortwave','latent','longwave','sensible') 
%       Qnet=latent+longwave+sensible+shortwave;
%         
%       % ##################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD climatology
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_clim;
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld = - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
%         Qnet      = Qnet      + shortwave_mld;
%         shortwave = shortwave + shortwave_mld;
%         clear shortwave_mld
%       % ##################################################################   
%        
% 
%       % ##################################################################
%       % Mixed layer warming due to Qnet
%         disp(['   Projected ML Warming in Year#',num2str(year)])
%         % Heat flux term in heat budget
%         MLT_Qnet_mon(:,:,:)  = Qnet     ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear Qnet
%         MLT_Qswr_mon(:,:,:)  = shortwave./mld_clim./3992./1027.*second_2_month; % K/month
%         clear shortwave
%         MLT_Qlat_mon(:,:,:)  = latent   ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear latent
%         MLT_Qlon_mon(:,:,:)  = longwave ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear longwave
%         MLT_Qsen_mon(:,:,:)  = sensible ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear sensible
%         % figure;imagesc(MLT_Qnet_mon(240:360,91:150,1)); caxis([-3 3])
% 
%         
%          % ################################################################
%          % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
%          dMLT_Qnet_mon0=MLT_Qnet_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qswr_mon0=MLT_Qswr_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlat_mon0=MLT_Qlat_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlon_mon0=MLT_Qlon_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qsen_mon0=MLT_Qsen_mon(240:360,91:150,:).*Sxy_NA(:,:);
% 
%          Sxy_NA0=Sxy_NA;
%          for month=2:12
%             Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
%          end
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          dMLT_Qnet_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlat_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlon_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qsen_mon0(isnan(Sxy_NA0))=NaN; 
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
% 
%          dMLT_Qnet_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qswr_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlat_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlat_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlon_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlon_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qsen_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qsen_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          clear dMLT_*_mon0 Sxy_NA0
%       % ###################################################################
%          
%          save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'_3_Qsw_prime.mat'],...
%               'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
%               'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%           clear NA_SST_mon* dMLT_*_mon
%          % ################################################################
%       count_yr=count_yr+1;
% end
% clear *_dail_clim
% % #########################################################################
% % #########################################################################  
% 
% 
% % #########################################################################
% % #########################################################################  
% %% 2.4 MLT Budget in 1981-2023 - Qsw_H' Term
% clc;clear
% time_ann=(1981:2010)';
% 
% % #########################################################################
%   % ACCESS Om2 0.25
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon lat
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% 
% % #########################################################################
%       % second_2_month=60*60*24*30;
%       second_2_month(1,1)  = 60*60*24*31;
%       second_2_month(2,1)  = 60*60*24*28.25;
%       second_2_month(3,1)  = 60*60*24*31;
%       second_2_month(4,1)  = 60*60*24*30;
%       second_2_month(5,1)  = 60*60*24*31;
%       second_2_month(6,1)  = 60*60*24*30;
%       second_2_month(7,1)  = 60*60*24*31;
%       second_2_month(8,1)  = 60*60*24*31;
%       second_2_month(9,1)  = 60*60*24*30;
%       second_2_month(10,1) = 60*60*24*31;
%       second_2_month(11,1) = 60*60*24*30;
%       second_2_month(12,1) = 60*60*24*31;
%       
%       for month=1:12
%           second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
%       end
%       second_2_month=second_2_month0; clear second_2_month0
% % #########################################################################
%       
% 
% % #########################################################################
%   % Monthly MLD Climatology from IAP
%     count=1;
%     for year=1981:2010
%         disp(['   MLD from IAPv3, Year#',num2str(year)])
%           % dRh0=0.125 ########################################################
%           load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%           lat_IAP=lat_IAP(:,1);
%           mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%           clear mld
% 
%           lon_IAP0(1:340,1)  = lon_IAP(21:360,1);
%           lon_IAP0(341:360,1)= lon_IAP(1:20,1)+360; clear lon_IAP
%           lon_IAP            = lon_IAP0;            clear lon_IAP0
% 
%           mld00(1:340,:,:)   = mld0(21:360,:,:);
%           mld00(341:360,:,:) = mld0(1:20,:,:); clear mld0
%           mld_monthly(:,:,:,count)=mld00;      clear mld00
%           count=count+1;
%     end
%     clear year count
%     mld_clim=nanmean(mld_monthly,4); clear mld_monthly
% % #########################################################################
% 
% 
% % #########################################################################
%   % Monthly SWR Climatology
%       count=1;
%       for year_qnet=1981:2010
%           disp(['    heat fluxes Year#',num2str(year_qnet)])
%           load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year_qnet),'_1deg.mat'],...
%                 'shortwave') 
%           shortwave0(:,:,:,count)=shortwave;
%           clear shortwave
%           count=count+1;
%       end
%       clear year_qnet count
%       SWR_clim=nanmean(shortwave0,4); clear shortwave0 
% % #########################################################################
%       
% 
% % #########################################################################
% count_yr=1;
% for year=1981:2023
%     disp([' Predicting NA daily SST Year#',num2str(year)])
%     % #####################################################################
%     % Monthly MLD from IAP
%       % dRh0=0.125 ####################################################
%       load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld')
%       mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%       clear mld
%       mld00(1:340,:,:)   = mld0(21:360,:,:);
%       mld00(341:360,:,:) = mld0(1:20,:,:); clear mld0
%       mld_mon(:,:,:)     = mld00;          clear mld00
%           
%           
%     % ##################################################################### 
%     % Monthly heat flux terms
%       load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
%             'shortwave','latent','longwave','sensible') 
%       Qnet=latent+longwave+sensible+SWR_clim;
%       
%       % ##################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_mon;
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld = - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
%         Qnet      = Qnet     + shortwave_mld;
%         shortwave = SWR_clim + shortwave_mld;
%         clear shortwave_mld
%       % ##################################################################   
%        
% 
%       % ##################################################################
%       % Mixed layer warming due to Qnet
%         disp(['   Projected ML Warming in Year#',num2str(year)])
%         % Heat flux term in heat budget
%         MLT_Qnet_mon(:,:,:)  = Qnet     ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear Qnet
%         MLT_Qswr_mon(:,:,:)  = shortwave./mld_clim./3992./1027.*second_2_month; % K/month
%         clear shortwave
%         MLT_Qlat_mon(:,:,:)  = latent   ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear latent
%         MLT_Qlon_mon(:,:,:)  = longwave ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear longwave
%         MLT_Qsen_mon(:,:,:)  = sensible ./mld_clim./3992./1027.*second_2_month; % K/month
%         clear sensible
%         % figure;imagesc(MLT_Qnet_mon(240:360,91:150,1)); caxis([-3 3])
% 
%         
%          % ################################################################
%          % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
%          dMLT_Qnet_mon0=MLT_Qnet_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qswr_mon0=MLT_Qswr_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlat_mon0=MLT_Qlat_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlon_mon0=MLT_Qlon_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qsen_mon0=MLT_Qsen_mon(240:360,91:150,:).*Sxy_NA(:,:);
% 
%          Sxy_NA0=Sxy_NA;
%          for month=2:12
%             Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
%          end
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          dMLT_Qnet_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlat_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlon_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qsen_mon0(isnan(Sxy_NA0))=NaN; 
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
% 
%          dMLT_Qnet_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qswr_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlat_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlat_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlon_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlon_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qsen_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qsen_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          clear dMLT_*_mon0 Sxy_NA0
%       % ###################################################################
%          
%          save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_4_Qsw_H_prime.mat'],...
%               'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
%               'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%           clear NA_SST_mon* dMLT_*_mon
%          % ################################################################
%       count_yr=count_yr+1;
% end
% clear *_dail_clim
% % #########################################################################
% % #########################################################################  




%% ########################################################################
%% Plotting 5a: Bar Charts for MLT Budget Decomposition 
clc;clear
time_ann=(1981:2022)';

% #########################################################################
% #########################################################################
% 1. MLTa from IAP data
% #########################################################################
%   % ACCESS OM2 0.25
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; 
       Sxy_NA=Sxy; clear Sxy % figure;imagesc(basin_mask_NA)
% #########################################################################



% #########################################################################
    % Monthly MLT during 1981-2023
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,1:12);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear month 
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count Sxy_NA lon lat

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(2:12,:)-mlt_month_19812023_ano(1:11,:);
    clear mlt_month* 
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.1 MLT Budget in 1981-2023
count_yr=1;
for year=1981:2023
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'.mat'],...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.2 MLT Budget in 1981-2023 - MLD' Term
count_yr=1;
for year=1981:2023
    % Test with MLD and Qsw_H climatology
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'_2_MLD_prime.mat'],...
          'dMLT_Qswr_mon')
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qswr_mon0,2)
        NA_MLT_Qswr_ano_2_MLD_prime(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.3 MLT Budget in 1981-2023 - Qsw' Term
count_yr=1;
for year=1981:2023
    % Test with MLD and Qsw_H climatology
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V4_3way_decomposition_',num2str(year),'_3_Qsw_prime.mat'],...
          'dMLT_Qswr_mon')
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qswr_mon0,2)
        NA_MLT_Qswr_ano_3_Qsw_prime(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.4 MLT Budget in 1981-2023 - Qsw_H' Term
count_yr=1;
for year=1981:2023
    % Test with MLD and SWR climatology
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_4_Qsw_H_prime.mat'],...
          'dMLT_Qswr_mon')
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qswr_mon0,2)
        NA_MLT_Qswr_ano_4_Qsw_H_prime(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.230; ixe = 0.230;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

pos{101}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 

clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);


subplot('position',pos{101}) 
   h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.94); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

       
       % ##################################################################
        % Bars and lines from MLD/Qsw_H/SWR climatology
        hold on
        bar_dmlt_Clim(1:9,1:21)=NaN;
        % MLD' Term
        bar_dmlt_Clim(1:9,7)=NA_MLT_Qswr_ano_2_MLD_prime(1:9,43);       
        % Qsw' Term
        bar_dmlt_Clim(1:9,8)=NA_MLT_Qswr_ano_3_Qsw_prime(1:9,43);
        % Qsw_H' Term
        bar_dmlt_Clim(1:9,9)=NA_MLT_Qswr_ano_4_Qsw_H_prime(1:9,43);  
        
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.68); 
           hold on
           set(h0_Clim(7),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(8),'FaceColor',[0.2, 0.3, 0.2],'EdgeColor',[0.2, 0.3, 0.2])
           hold on
           set(h0_Clim(9),'FaceColor',[0.2, 0.8, 0.6],'EdgeColor',[0.2, 0.8, 0.6])
           

           hold on
           plot((6.10:0.001:6.23),0.70*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.10:0.001:6.23),0.59*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.10:0.001:6.23),0.48*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)  
           hold on
           plot((6.10:0.001:6.23),0.70*ones(length((6.10:0.001:6.23)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.10:0.001:6.23),0.59*ones(length((6.10:0.001:6.23)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.10:0.001:6.23),0.48*ones(length((6.10:0.001:6.23)')),'-','color',[0.2, 0.8, 0.6],'linewidth',4)  

           text(6.25, 0.70,'$\mathbf{MLD^\prime}$', 'Interpreter', 'latex','fontsize',18,'FontName', 'Aptos')
           text(6.25, 0.59,'$\mathbf{{Q_{sw}}^\prime}$', 'Interpreter', 'latex','fontsize',18,'FontName', 'Aptos')
           text(6.25, 0.48,'$\mathbf{{Q_{sw,H}}^\prime}$', 'Interpreter', 'latex','fontsize',18,'FontName', 'Aptos')
       % ##################################################################      
       
       
    % #####################################################################
    hold on
    leg1=legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2);
    hold on
    set(leg1,'fontsize',18)
    hold on
    legend('boxoff')
    % #####################################################################


    % Legend will show names for each color
    set(gca,'Ylim',[-0.44 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.6)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2','1.6'},'fontsize',21)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:1:8.5)
    set(gca,'XTickLabel',{[],[],[],[]},'fontsize',21)

    text(4.915,-0.54,'May','fontsize',21,'color','k','FontWeight','normal')
    text(5.930,-0.54,'Jun','fontsize',21,'color','k','FontWeight','normal')
    text(6.945,-0.54,'Jul','fontsize',21,'color','k','FontWeight','normal')
    text(7.920,-0.54,'Aug','fontsize',21,'color','k','FontWeight','normal')

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
    ylabel(['[ \circC per month ]'],'fontsize',21,'color','k','FontWeight','normal')
    
    text(4.5,0.88,['(a) Observed mixed layer temperature budget anomalies'],'fontsize',22,'color','k','FontWeight','bold')

% #########################################################################
% #########################################################################  

disp(['>> !!'])
disp(['  Save to left built-in screen...'])
disp(['>> !!'])


