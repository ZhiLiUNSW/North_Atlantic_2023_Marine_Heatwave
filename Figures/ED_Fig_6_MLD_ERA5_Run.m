%% Extended Data Figure: Monthly average model simulated MLD anomalies in 2023
%  May, June, July, August, Monthly Mean 1981-2023
%  Monthly MLD from ACCESS OM2, 1981-2023


%% ########################################################################
%  0.1 MLD anomaly from ERA5-forced ACCESS OM2 simulation
clc; clear
time_ann(:,1)=(1980:1:2023);    
time_mon=(1980:1/12:2011-1/12)';

% #########################################################################
      % Basin Mask from ACCESS OM2 0.25 forced by ERA5
      load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
      [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
      Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
      Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
% #########################################################################
% MLD for 1980-2010
    % #####################################################################
    % Monthly MLD from ACCESS, 1980-2010
    disp(['   MLD from ERA5-forced run, Year#1980-2010'])
     lon=ncread('MLD_025deg_era5_iaf_1958cycle1_1980-2010_clim.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_era5_iaf_1958cycle1_1980-2010_clim.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_era5_iaf_1958cycle1_1980-2010.nc','mld');
     mld=mld0(721:1200,499:891,:); clear mld0
     

     count_yr=1;
     for year=1980:2010
         disp(['   MLT budget in Year#',num2str(year)])
         mld_monthly=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);

         % ################################################################
         % NA (100W-20E,0-60N) averaged MLD
         mld_monthly=mld_monthly(:,1:302,:).*Sxy_NA(:,1:302);
         Sxy_NA0=Sxy_NA(:,1:302,1);
         for month=2:12
            Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
         end
         clear month
         Sxy_NA0(isnan(mld_monthly))=NaN; 
         mld_monthly(isnan(Sxy_NA0))=NaN; 

         mld_mon(1:12,count_yr)=squeeze(nansum(nansum(mld_monthly(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         clear dMLT_*_mon0 Sxy_NA0
         % ###################################################################
         count_yr=count_yr+1;
     end
     clear year mld_monthly count_yr mld
     
     mld_mon_19812010_clim=squeeze(nanmean(mld_mon(:,1:30),2));
% #########################################################################
% #########################################################################  
  
     
% #########################################################################
% MLD for 2023
    % #####################################################################
    % Monthly MLD from ACCESS, 2023
    disp(['   MLD from ERA5-forced run, Year#2023'])
     lon=ncread('MLD_025deg_era5_iaf_1958cycle1_2023.nc','xt_ocean'); 
     lon=lon(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_era5_iaf_1958cycle1_2023.nc','yt_ocean'); 
     lat=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_era5_iaf_1958cycle1_2023.nc','mld');
     mld(:,:,1:12)=mld0(721:1200,499:891,1:12); clear mld0
     
     % ####################################################################
     % NA (100W-20E,0-60N) averaged MLD
     mld_monthly=mld(:,1:302,:).*Sxy_NA(:,1:302);
     Sxy_NA0=Sxy_NA(:,1:302,1);
     for month=2:12
        Sxy_NA0(:,:,month)=Sxy_NA0(:,1:302,1); 
     end
     clear month
     Sxy_NA0(isnan(mld_monthly))=NaN; 
     mld_monthly(isnan(Sxy_NA0))=NaN; 

     mld_2023(1:12,1)=squeeze(nansum(nansum(mld_monthly(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
     clear dMLT_*_mon0 Sxy_NA0
% #########################################################################


% #########################################################################
%  Plotting: Normalized MLD anomaly, 1*STD
%   Anomaly is relative to the clim-daily mean
    count=1;
    for year=1:size(mld_mon,2)
        mld_monthly_ano(:,count)=mld_mon(:,count)-mld_mon_19812010_clim;
        count=count+1;
    end    
    mld_mon_2023_ano(:,1)=mld_2023(:,1)-mld_mon_19812010_clim;
    
    % Normalize the anomalies
    for month=1:12
        mld_monthly_std(month,1)=std(squeeze(mld_monthly_ano(month,1:30)));
    end
    
    mld_monthly_ano   = mld_monthly_ano  ./mld_monthly_std ./1;
    mld_mon_2023_ano  = mld_mon_2023_ano ./mld_monthly_std ./1;
    
    save('plot_2023_MLD_Anomaly_RG18_2_Spagatti_V7_1_ACCESS_ERA5_MLD.mat','mld_mon_2023_ano','mld_monthly_ano')
% #########################################################################
% #########################################################################



%% ########################################################################
%  0.2 MLD anomaly from JRA55-Forced ACCESS OM2 simulation    
clc; clear
time_ann(:,1)=(1981:1:2023);    
time_mon=(1981:1/12:2011-1/12)';

% #########################################################################
      % Basin Mask from ACCESS OM2 0.25 forced by ERA5
      load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
      [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
      Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
      Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
clear mld*
% Monthly MLD during 1981-2023, ACCESS OM2
    time_ann(:,1)=(1981:1:2023);    
    count=1;
    for year=1981:2023
          disp(['loading mld... year# ',num2str(year)])
          if year==2018
              mld00(:,:,1:12)=NaN;
          else
              mld00(:,:,1:12)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mld');
          end
          
          mld00=mld00(:,:,:).*Sxy_NA;
          Sxy_NA0=Sxy_NA;
          Sxy_NA0(isnan(mld00(:,:,1)))=NaN;

          % NA 0-60N
          mld_monthly(1:12,count)=squeeze(nansum(nansum(mld00(:,1:302,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:302),2),1));
          clear mld00 Sxy_NA0
          
          count=count+1;
    end
    clear year count
    mld_monthly(mld_monthly==0)=NaN;
% #########################################################################


% #########################################################################
%  Plotting: Normalized MLD anomaly, 1*STD
%  Anomaly is relative to the clim-daily mean
    mld_mon_19812010_clim=squeeze(nanmean(mld_monthly(:,1:30),2));
    count=1;
    for year=1:size(mld_monthly,2)
        mld_monthly_ano(:,count)=mld_monthly(:,count)-mld_mon_19812010_clim;
        count=count+1;
    end    
    mld_mon_2023_ano(:,1)=mld_monthly_ano(1:12,end);
    
    % Normalize the anomalies
    for month=1:12
        mld_monthly_std(month,1)=std(squeeze(mld_monthly(month,1:30)));
    end
    
    mld_monthly_ano   = mld_monthly_ano  ./mld_monthly_std ./1;
    mld_mon_2023_ano  = mld_mon_2023_ano ./mld_monthly_std ./1;
    

    save('plot_2023_MLD_Anomaly_RG18_2_Spagatti_V7_2_ACCESS_JRA55_MLD.mat','mld_mon_2023_ano','mld_monthly_ano')
% #########################################################################
% #########################################################################




%% ########################################################################
%% Mapping 2023 MLD Anomaly Relative to 1981-2010 Mean - ACCESS OM2, ERA5&JRA55
clc;clear
     

% #########################################################################
% #########################################################################
% MLD for 1980-2010
    % #####################################################################
    % Monthly MLD from ACCESS, 1980-2010
    disp(['   MLD from ERA5-forced run, Year#1980-2010'])
     lon_025=ncread('MLD_025deg_era5_iaf_1958cycle1_1980-2010_clim.nc','xt_ocean'); 
     lon_025=lon_025(721:1200,1)+360; % 100W-20E
     lat=ncread('MLD_025deg_era5_iaf_1958cycle1_1980-2010_clim.nc','yt_ocean'); 
     lat_025=lat(499:891,1);% 0-70N
     
     mld0=ncread('MLD_025deg_era5_iaf_1958cycle1_1980-2010.nc','mld');
     mld=mld0(721:1200,499:891,:); clear mld0

     count_yr=1;
     for year=1980:2010
         disp(['   MLT budget in Year#',num2str(year)])
         mld_monthly(:,:,1:12,count_yr)=mld(:,:,(count_yr-1)*12+1:(count_yr-1)*12+12);
         count_yr=count_yr+1;
     end
     clear year count_yr mld
% #########################################################################
% MLD for 2023
    % #####################################################################
    % Monthly MLD from ACCESS, 2023
    disp(['   MLD from ERA5-forced run, Year#2023'])
     mld0=ncread('MLD_025deg_era5_iaf_1958cycle1_2023.nc','mld');
     mld_mon_2023(:,:,1:12)=mld0(721:1200,499:891,1:12); clear mld0
 
    % Anomaly is relative to the clim-daily mean
    mld_mon_19812010_clim=squeeze(nanmean(mld_monthly,4));
    mld_mon_2023_ano(:,:,:)=mld_mon_2023-mld_mon_19812010_clim;

    mld_mon_19812010=squeeze(mld_monthly(:,:,:,1:31));
    % Normalize the anomalies
    for i=1:size(mld_mon_19812010,1)
        for j=1:size(mld_mon_19812010,2)
            for month=1:12
                mld_mon_19812010_std(i,j,month,1)=std(squeeze(mld_mon_19812010(i,j,month,1:end)));
            end
        end
    end

    mld_mon_2023_ano=mld_mon_2023_ano./mld_mon_19812010_std;
% #########################################################################
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
% #########################################################################

    
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,smooth2a(mld_mon_2023_ano(:,:,5)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('a. MLD anomaly (May 2023)','fontsize',20,'FontWeight','bold')
        
 

subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,smooth2a(mld_mon_2023_ano(:,:,6)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('b. MLD anomaly (June 2023)','fontsize',20,'FontWeight','bold')
        
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,smooth2a(mld_mon_2023_ano(:,:,7)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('c. MLD anomaly (July 2023)','fontsize',20,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,smooth2a(mld_mon_2023_ano(:,:,8)',0,0));
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        %caxis([-24 24]);
        title('d. MLD anomaly (August 2023)','fontsize',20,'FontWeight','bold')
        

        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.012 iys+1.5*iyw+0*iyd 0.010 1*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.05);
        ylabel(hBar1, '[ 1 std. dev. ]','rotation',90);
        
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
%  5. Time series of normalized MLD from ACCESS OM2 ERA5 and JRA55
clc;
clear mld* basin*
pos{51}  = [ixs+0.0*ixw  iys+0*iyw+0*iyd   2.0*ixw+1*ixd 0.918*iyw]; 

% #########################################################################
% MLD anoamly from ACCESS OM2 ERA5
    load('plot_2023_MLD_Anomaly_RG18_2_Spagatti_V7_1_ACCESS_ERA5_MLD.mat','mld_mon_2023_ano','mld_monthly_ano')
    mld_2023_access_era5 = mld_mon_2023_ano; clear mld_mon_2023_ano
     mld_mon_access_era5 = mld_monthly_ano;  clear mld_monthly_ano
% #########################################################################


% #########################################################################
% MLD anoamly from ACCESS OM2 JRA55
    load('plot_2023_MLD_Anomaly_RG18_2_Spagatti_V7_2_ACCESS_JRA55_MLD.mat','mld_mon_2023_ano','mld_monthly_ano')
    mld_2023_access_jra55 = mld_mon_2023_ano; clear mld_mon_2023_ano
     mld_mon_access_jra55 = mld_monthly_ano;  clear mld_monthly_ano
% #########################################################################
    


% Plotting ################################################################
clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
subplot('position',pos{51})
    for year=1:size(mld_mon_access_era5,2)
        line00=plot(1:12,runmean(squeeze(mld_mon_access_era5(:,year)),0));
        set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
        hold on
    end
    line00=plot(1:12,runmean(squeeze(mld_mon_access_era5(:,5)),0));
    set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 


% #########################################################################
    % 2023 MLD anomaly from ACCESS OM2 JRA55
    hold on
    line_jra55=plot(1:12,runmean(squeeze(mld_2023_access_jra55(:,1)),0));
    set(line_jra55,'color',[0.99, 0.794, 0.025],'LineWidth',5,'linestyle','-'); 
    % 2023 MLD anomaly from ACCESS OM2 ERA5
    hold on
    line_era5=plot(1:12,runmean(squeeze(mld_2023_access_era5(:,1)),0));
    set(line_era5,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
% #########################################################################

    
    leg=legend([line00 line_era5 line_jra55],'1981-2010','2023, ACCESS OM2, ERA5','2023, ACCESS OM2, JRA55',...
               'Location','northwest','NumColumns',2);
    set(leg,'fontsize',18)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
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
    ylabel(['[ 1 std. dev. ]'],'fontsize',20,'color','k','FontWeight','normal')
    % xlabel(['Month'],'fontsize',24,'color','k','FontWeight','normal')
    
    title('e. Normalized MLD anomaly in 1981-2023, ACCESS OM2','fontsize',20,'color','k','FontWeight','bold')

% #########################################################################
% #########################################################################
   
disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')



