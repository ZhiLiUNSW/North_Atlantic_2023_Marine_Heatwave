%% Main Figure 2: NA MLD Anomaly in 2023 Relative to 1981-2010 Mean
%  May, June, July, August, in 2023
%  Monthly MLD from IAP TS, 1981-2023


%% #######################################################################
%% Figure 1# 2023 MLD Anomaly Relative to 1981-2010 Mean - IAP data
clc;clear

% #########################################################################
% 1. Monthly MLD during 1981-2010 - from IAP
count=1;
time_ann(:,1)=(1981:1:2023);
for year=1981:2023
      disp(['loading mld... year# ',num2str(year)])
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',                  num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,:)=squeeze(mld(:,:,1,1:12)); % 64.5S-64.5N
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0

      mld_monthly(:,:,1:12,count)=squeeze(mld00);
      clear mld00 
      count=count+1;
end
clear year count
mld_monthly(mld_monthly==0)=NaN;

    % Anomaly is relative to the clim-daily mean
    mld_mon_19812010_clim=squeeze(nanmean(mld_monthly(:,:,:,1:30),4));
    mld_mon_2023_ano(:,:,:)=mld_monthly(:,:,:,end)-mld_mon_19812010_clim;
    
    % Normalize the anomalies
    mld_mon_19812010=squeeze(mld_monthly(:,:,:,1:30));
    for i=1:size(mld_mon_19812010,1)
        for j=1:size(mld_mon_19812010,2)
            for month=1:12
                mld_mon_19812010_std(i,j,month,1)=std(squeeze(mld_mon_19812010(i,j,month,1:30)));
            end
        end
    end

    mld_mon_2023_ano=mld_mon_2023_ano./mld_mon_19812010_std;
% #########################################################################



% #########################################################################
% 2. SST Anomaly Relative Clim Mean State
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

           

% #########################################################################
% Plotting ################################################################
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.270; ixe = 0.270;  ixd = 0.030; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.040; iye = 0.040;  iyd = 0.0450; iyw = (1-iys-iye-2*iyd)/3;

pos{11}  = [ixs          iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=color(14:-1:9,:);
color0(7:12,:)=color(7:-1:2,:);

    
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_mon_2023_ano(:,:,5)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.25);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_May,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.25);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('a. MLD anomaly (May 2023)','fontsize',19,'FontWeight','bold')
        
 

subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_mon_2023_ano(:,:,6)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.25);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jun,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.25);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('b. MLD anomaly (June 2023)','fontsize',19,'FontWeight','bold')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_mon_2023_ano(:,:,7)',0,0));
        shading flat
        hold on

        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.25);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Jul,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.25);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('c. MLD anomaly (July 2023)','fontsize',19,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(mld_mon_2023_ano(:,:,8)',0,0));
        shading flat
        hold on
        
        % #########################################################################
        % SST anomaly
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[-1 -1],'color',[0.301, 0.745, 0.933],'linewidth',3.25);
        v21=[];
        clabel(C21,h21,v21,'labelspacing',18000,'fontsize',14)
        hold on;
        [C21,h21]=m_contour(lon_ERA5,lat_ERA5,smooth2a(SST_ano_Aug,12,12)',[1 1],'color',[0.850, 0.325, 0.098],'linewidth',3.25);
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
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('d. MLD anomaly (August 2023)','fontsize',19,'FontWeight','bold')
        

        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.012 iys+1.5*iyw+0*iyd 0.010 1*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',19,'FontName','Arial','LineWidth',1.2,'TickLength',0.05);
        ylabel(hBar1, '[ 1 std. dev. ]','rotation',90);
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% Twitter part
clc;
clear mld* basin*
pos{51}  = [ixs+0.0*ixw  iys+0*iyw+0*iyd   2.0*ixw+1*ixd 0.918*iyw]; 

% #########################################################################
% MLD anoamly from ACCESS OM2
    load('plot_2023_MLD_Anomaly_RG18_2_Spagatti_V6_IAP_and_ACCESS_MLD_3_ACCESS.mat','mld_mon_2023_ano')
    mld_2023_access=mld_mon_2023_ano; clear mld_mon_2023_ano
% #########################################################################

% #########################################################################
% Basin Mask from ACCESS OM2 0.25
     clear basin* Sxy*
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
% Monthly MLD during 1981-2010
    count=1;
    time_ann(:,1)=(1981:1:2023);
    for year=1981:2023
          disp(['loading mld... year# ',num2str(year)])
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(26:155,1);
          if year==2023
              mld0(:,:,:,1:9)=mld(:,:,:,1:9); clear mld
              load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAPv4_Temp_2023.mat'],'mld','lon_IAP','lat_IAP')
              mld0(:,:,:,10:12)=mld(:,:,:,10:12); clear mld
              mld=mld0; clear mld0
          end
          mld0(:,:,:)=squeeze(mld(:,26:155,1,1:12)); % 64.5S-64.5N
          clear mld
          
          lon_IAP0(1:340,1)=lon_IAP(21:360,1);
          lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
          lon_IAP=lon_IAP0; clear lon_IAP0

          mld00(1:340,:,:)=mld0(21:360,:,:);
          mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0

          mld00=mld00(240:360,66:125,:).*Sxy_NA;
          Sxy_NA0=Sxy_NA;
          Sxy_NA0(isnan(mld00(:,:,1)))=NaN;
          lon_IAP0=lon_IAP(240:360,1);

          % NA 0-60N
          mld_monthly(1:12,count)=squeeze(nansum(nansum(mld00(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,:),2),1));
          % Eastern and Western NA 0-60N, 40W
          mld_monthly_ENA(1:12,count)=squeeze(nansum(nansum(mld00(62:end,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(62:end,:),2),1));
          mld_monthly_WNA(1:12,count)=squeeze(nansum(nansum(mld00(1:61,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(1:61,:),2),1));
          clear mld00 Sxy_NA0
          count=count+1;
    end
    clear year count
      mld_monthly(mld_monthly==0)=NaN;
      mld_monthly_ENA(mld_monthly_ENA==0)=NaN;
      mld_monthly_WNA(mld_monthly_WNA==0)=NaN;
% #########################################################################


% #########################################################################
%  Plotting: Normalized MLD anomaly, 1*STD
%   Anomaly is relative to the clim-daily mean
    mld_mon_19812010_clim=squeeze(nanmean(mld_monthly(:,1:30),2));
    mld_mon_19812010_clim_ENA=squeeze(nanmean(mld_monthly_ENA(:,1:30),2));
    mld_mon_19812010_clim_WNA=squeeze(nanmean(mld_monthly_WNA(:,1:30),2));
    count=1;
    for year=1:size(mld_monthly,2)
        mld_monthly_ano(:,count)=mld_monthly(:,count)-mld_mon_19812010_clim;
        mld_monthly_ano_ENA(:,count)=mld_monthly_ENA(:,count)-mld_mon_19812010_clim_ENA;
        mld_monthly_ano_WNA(:,count)=mld_monthly_WNA(:,count)-mld_mon_19812010_clim_WNA;
        count=count+1;
    end    
    mld_mon_2023_ano(:,1)=mld_monthly_ano(1:12,end);
    mld_mon_2023_ano_ENA(:,1)=mld_monthly_ano_ENA(1:12,end);
    mld_mon_2023_ano_WNA(:,1)=mld_monthly_ano_WNA(1:12,end);

    
    % Normalize the anomalies
    for month=1:12
        mld_monthly_std(month,1)=std(squeeze(mld_monthly(month,1:30)));
        mld_monthly_std_ENA(month,1)=std(squeeze(mld_monthly_ENA(month,1:30)));
        mld_monthly_std_WNA(month,1)=std(squeeze(mld_monthly_WNA(month,1:30)));
    end
    
    mld_monthly_ano=mld_monthly_ano./mld_monthly_std./1;
    mld_mon_2023_ano    =mld_mon_2023_ano    ./mld_monthly_std./1;
    mld_mon_2023_ano_ENA=mld_mon_2023_ano_ENA./mld_monthly_std_ENA./1;
    mld_mon_2023_ano_WNA=mld_mon_2023_ano_WNA./mld_monthly_std_WNA./1;


% Plotting ################################################################
clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
subplot('position',pos{51})
    for year=1:42
        line00=plot(1:12,mld_monthly_ano(:,year));
        set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
        hold on
    end
    line00=plot(1:12,mld_monthly_ano(:,5));
    set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 


% #########################################################################
    % 2023 MLD anomaly from ACCESS OM2
    % 2023 MLD anomaly from IAP
    hold on
    line23=plot(1:12,mld_mon_2023_ano(:,1));
    set(line23,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
    hold on
    line_WNA=plot(1:12,mld_mon_2023_ano_WNA(:,1));
    set(line_WNA,'color',[0.950,0.325,0.098],'LineWidth',2.5,'linestyle','--'); 
    hold on
    line_ENA=plot(1:12,mld_mon_2023_ano_ENA(:,1));
    set(line_ENA,'color',[0.950,0.325,0.098],'LineWidth',2.5,'linestyle',':'); 
% #########################################################################

    
    leg=legend([line00 line23 line_ENA line_WNA],'1981-2022','2023','2023, Eastern NA','2023, Western NA','Location','northwest','NumColumns',2);
    set(leg,'fontsize',18)
    hold on
    legend('boxoff')

    
    set(gca,'Ylim',[-4 5],'ycolor','k') 
    set(gca,'YTick',-6:2:6)
    set(gca,'YTickLabel',{'-6','-4','-2','0','2','4','6'},'fontsize',19)
    set(gca,'Xlim',[1 12]) 
    set(gca,'XTick',1:1:12)
    set(gca,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun',...
                          'Jul','Aug','Sep','Oct','Nov','Dec'},'fontsize',19)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ 1 std. dev. ]'],'fontsize',19,'color','k','FontWeight','normal')
    
    title('e. Normalized MLD Anomaly in 1981-2023','fontsize',19,'color','k','FontWeight','bold')

% ########################################################################
disp(' Save to the built-in screen')
disp(' >>')

