%% Extended Data Figure: Trends in summertime MLD, SST, SSS and SSD during 1981-2023
%  June, July, August Averaged
%  Monthly Fields from IAP SSS and MLD, and ERA5 SST, 1981-2023


%% Plotting SST, SSS, MLD trend ###########################################
%  June, July, August Averaged
clc;clear
time_ann(:,1)=(1981:1:2023);


% #########################################################################
% 1. June-August SSS during 1981-2023 - from IAP
% Monthly SSS during 1981-2023
count=1;
time_ann(:,1)=(1981:1:2023);
for year=1981:2023
      disp(['loading SSS... year# ',num2str(year)])
      load(['SA_depth_IAP_C17_',num2str(year),'.mat'],'SA','lon_IAP','lat_IAP')
      lat_IAP=lat_IAP(:,1);
      SA_depth0(:,:,:)=squeeze(nanmean(SA(1,:,:,1:12),1));
      clear SA

      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0

      SA_depth00(1:340,:,:)=SA_depth0(21:360,:,:);
      SA_depth00(341:360,:,:)=SA_depth0(1:20,:,:); clear SA_depth0

      SSS_mon_19812023(:,:,1:12,count)=squeeze(SA_depth00);
      clear SA_depth00
      count=count+1;
end
clear year count
SSS_mon_19812023(SSS_mon_19812023==0)=NaN;


SSS_ann_19812023=squeeze(nanmean(SSS_mon_19812023(:,:,6:8,:),3)); % Jun-Aug Averaged
clear SSS_mon_19812023

% Spatial Trend 1981-2023
for i=1:size(SSS_ann_19812023,1)
    disp(['  spatial trend#SSS lon#',num2str(i)])
    for j=1:size(SSS_ann_19812023,2)
        p1=polyfit(1:size(SSS_ann_19812023,3),squeeze(SSS_ann_19812023(i,j,1:end))',1);
        SSS_trend(i,j)=p1(1); 
        clear p1
    end
end
clear i j
% #########################################################################


% #########################################################################
% 2. June-August MLD during 1981-2023 - from IAP
count=1;
for year=1981:2023
      disp(['  loading mld... year# ',num2str(year)])
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
      if size(mld,4)<12
          mld(:,:,:,size(mld,4)+1:12)=NaN;
      end
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,:)=squeeze(mld(:,:,1,1:12)); % 64.5S-64.5N
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      % mld00(isnan(basin_dep_mon))=NaN;
      
      mld_mon_19812023(:,:,1:12,count)=squeeze(mld00);
      clear mld00 
      count=count+1;
end
clear year count
mld_mon_19812023(mld_mon_19812023==0)=NaN;
mld_mon=squeeze(nanmean(mld_mon_19812023(:,:,6:8,:),3));
clear mld_mon_19812023

% Spatial Trend 1981-2023
for i=1:size(mld_mon,1)
    disp(['  spatial trend#MLD lon#',num2str(i)])
    for j=1:size(mld_mon,2)
        p1=polyfit(1:size(mld_mon,3),squeeze(mld_mon(i,j,1:end))',1);
        MLD_trend(i,j)=p1(1); 
        clear p1
    end
end
clear i j
% #########################################################################


% #########################################################################
% 3. June-August SST during 1981-2023 - from ERA5
count=1;
for year=1981:2023
      disp(['  loading SST... year# ',num2str(year)])
      load(['SST_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'sst','lon','lat')
      SST0(:,:,:,count)=squeeze(sst(:,1:end,:));
      clear sst
      count=count+1;
end
clear year count
SST_mon=squeeze(nanmean(SST0(:,:,6:8,:),3));
clear SST0

lon0(1:1360,1)=lon(81:1440,1);
lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
lon_ERA5=lon0; clear lon0
lat_ERA5=lat;  clear lat

SST_mon0(1:1360,:,:)=SST_mon(81:1440,:,:);
SST_mon0(1361:1440,:,:)=SST_mon(1:80,:,:); clear SST_mon
SST_mon=SST_mon0; clear SST_mon0

% Spatial Trend 1981-2023
for i=1:size(SST_mon,1)
    disp(['  spatial trend#SST lon#',num2str(i)])
    for j=1:size(SST_mon,2)
        p1=polyfit(1:size(SST_mon,3),squeeze(SST_mon(i,j,1:end))',1);
        SST_trend(i,j)=p1(1); 
        clear p1
    end
end
clear i j
% #########################################################################
 


% #########################################################################
% 4.1 June, July, August Sfc Density during 1981-2023 - from IAP
load(['plot_2023_Sfc_Density_by_CT_SA_IAP_Monthly_V1_1_RealT_RealS.mat'],'Rh0_depth','lon_IAP','lat_IAP','depth10','time_ann');
Den_mon=squeeze(nanmean(Rh0_depth(:,:,1,6:8,:),4)); % Jun-Aug Averaged
clear Rh0_depth

lon_IAp0(1:340,1)=lon_IAP(21:360,1);
lon_IAp0(341:360,1)=lon_IAP(1:20,1)+360;
lon_IAP=lon_IAp0; clear lon_IAp0

Den_mon0(1:340,:,:)=Den_mon(21:360,:,:);
Den_mon0(341:360,:,:)=Den_mon(1:20,:,:);
Den_mon=Den_mon0; clear Den_mon0

% Spatial Trend 1981-2023
for i=1:size(Den_mon,1)
    disp(['  spatial trend#Dens lon#',num2str(i)])
    for j=1:size(Den_mon,2)
        p1=polyfit(1:size(Den_mon,3),squeeze(Den_mon(i,j,1:end))',1);
        Dens_trend(i,j)=p1(1); 
        clear p1
    end
end
clear i j Den_mon
% #########################################################################


% #########################################################################
% 4.2. June-August Sfc Density by CT during 1981-2023 - from IAP
load(['plot_2023_Sfc_Density_by_CT_SA_IAP_Monthly_V1_2_RealT_ClimS.mat'],'Rh0_depth');
Den_mon=squeeze(nanmean(Rh0_depth(:,:,1,6:8,:),4)); % Jun-Aug Averaged
clear Rh0_depth

Den_mon0(1:340,:,:)=Den_mon(21:360,:,:);
Den_mon0(341:360,:,:)=Den_mon(1:20,:,:);
Den_mon=Den_mon0; clear Den_mon0

% Spatial Trend 1981-2023
for i=1:size(Den_mon,1)
    disp(['spatial trend#Dens-CT lon#',num2str(i)])
    for j=1:size(Den_mon,2)
        p1=polyfit(1:size(Den_mon,3),squeeze(Den_mon(i,j,1:end))',1);
        Dens_trend_CT(i,j)=p1(1); 
        clear p1
    end
end
clear i j Den_mon
% #########################################################################


% #########################################################################
% 4.3. June-August Sfc Density by SA during 1981-2023 - from IAP
load(['plot_2023_Sfc_Density_by_CT_SA_IAP_Monthly_V1_3_ClimT_RealS.mat'],'Rh0_depth');
Den_mon=squeeze(nanmean(Rh0_depth(:,:,1,6:8,:),4)); % Jun-Aug Averaged
clear Rh0_depth

Den_mon0(1:340,:,:)=Den_mon(21:360,:,:);
Den_mon0(341:360,:,:)=Den_mon(1:20,:,:);
Den_mon=Den_mon0; clear Den_mon0

% Spatial Trend 1981-2023
for i=1:size(Den_mon,1)
    disp(['spatial trend#Dens-SA lon#',num2str(i)])
    for j=1:size(Den_mon,2)
        p1=polyfit(1:size(Den_mon,3),squeeze(Den_mon(i,j,1:end))',1);
        Dens_trend_SA(i,j)=p1(1); 
        clear p1
    end
end
clear i j Den_mon
% #########################################################################



% #########################################################################
%% Trend map for June-August SST SSS MLD
figure('Color',[1 1 1]);  %create a new figure of white color background

ixs = 0.220; ixe = 0.220;  ixd = 0.040; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.440; iye = 0.030;  iyd = 0.120; iyw = (1-iys-iye-1*iyd)/2;

pos{11}  = [ixs+0*ixw+0*ixd   iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+1*ixw+1*ixd   iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs+2*ixw+2*ixd   iys+1*iyw+1*iyd   ixw 1.0*iyw]; 

pos{12}  = [ixs+0*ixw+0*ixd   iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd   iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{32}  = [ixs+2*ixw+2*ixd   iys+0*iyw+0*iyd   ixw 1.0*iyw]; 

clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=color(14:-1:9,:);
color0(7:12,:)=color(7:-1:2,:);

% #########################################################################
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(MLD_trend(:,:)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.6 0.6]);
        title('a. Summertime MLD trend','fontsize',16,'FontWeight','bold')
        
        
        hBar3 = colorbar('EastOutside','horizontal');

        get(hBar3, 'Position');
        set(hBar3, 'Position', [ixs+0.05*ixw+0*ixd iys+1*iyw+0.82*iyd 0.9*ixw+0*iyd 0.010]);
        set(hBar3, 'ytick',-0.6:0.1:0.6,'yticklabel',{'<-0.6',[],'-0.4',[],'-0.2',[],'0',[],'0.2',[],'0.4',[],'>0.6'},'fontsize',15,'FontName','Arial','LineWidth',1.2,'TickLength',0.040);
        ylabel(hBar3, '[ m year^{-1} ]','rotation',0);   
    
    

subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_ERA5,lat_ERA5,smooth2a(SST_trend(:,:)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.06 0.06]);
        title('b. Summertime SST trend','fontsize',16,'FontWeight','bold')
        
         
        hBar1 = colorbar('EastOutside','horizontal');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+1.05*ixw+1*ixd iys+1*iyw+0.82*iyd 0.9*ixw+0*iyd 0.010]);
        set(hBar1, 'ytick',-0.06:0.01:0.06,'yticklabel',{'<-0.06',[],[],'-0.03',[],[],'0',[],[],'0.03',[],[],'>0.06'},'fontsize',15,'FontName','Arial','LineWidth',1.2,'TickLength',0.040);
        ylabel(hBar1, '[ \circC year^{-1} ]','rotation',0);   
    


subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(SSS_trend(:,:)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.012 0.012]);
        title('c. Summertime SSS trend','fontsize',16,'FontWeight','bold')
        
         
        hBar2 = colorbar('EastOutside','horizontal');

        get(hBar2, 'Position');
        set(hBar2, 'Position', [ixs+2.05*ixw+2*ixd iys+1*iyw+0.82*iyd 0.9*ixw+0*iyd 0.010]);
        set(hBar2, 'ytick',-0.012:0.002:0.012,'yticklabel',{'<-0.012',[],[],'-0.006',[],[],'0',[],[],'0.006',[],[],'>0.012'},'fontsize',15,'FontName','Arial','LineWidth',1.2,'TickLength',0.040);
        ylabel(hBar2, '[ g kg^{-1} year^{-1} ]','rotation',0);   
% #########################################################################



% #########################################################################
% Trend map for June-Auguts Sfc Density by CT and SA       
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(Dens_trend(:,:)',0,0));
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.012 0.012]);
        title('d. Summer surface density trend','fontsize',16,'FontWeight','bold')
        
         

subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(Dens_trend_CT(:,:)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.012 0.012]);
        title('e. Density trend due to temperature','fontsize',16,'FontWeight','bold')
        
        hBar4 = colorbar('EastOutside','horizontal');

        get(hBar4, 'Position');
        set(hBar4, 'Position', [ixs+1*ixw+1*ixd iys+0*iyw-0.18*iyd 1*ixw+0*iyd 0.010]);
        set(hBar4, 'ytick',-0.012:0.002:0.012,'yticklabel',{'<-0.012',[],[],'-0.006',[],[],'0',[],[],'0.006',[],[],'>0.012'},'fontsize',15,'FontName','Arial','LineWidth',1.2,'TickLength',0.030);
        %ylabel(hBar4, '[kg m^{-3} year^{-1}]','rotation',0);   
        m_text(395,-7,'[ kg m^{-3} year^{-1} ]','fontsize',16); 
     
        
        
subplot('position',pos{32})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(Dens_trend_SA(:,:)',0,0));
        shading flat
        hold on

        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',16,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.40 0.40 0.40]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.012 0.012]);
        title('f. Density trend due to salinity','fontsize',16,'FontWeight','bold')
        
% #########################################################################
        


% #########################################################################
% Vertical Profiles of NA CT, in (a-d) May-Aug
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/1_SST_Temp')
    load('plot_2023_SST_Anomaly_RG18_3_Temp_Profiles_Spagatti_V3_from19812010_IAP.mat','CT_ano_NA','CT_ano_ENA','CT_ano_WNA')
    clear depth10
    depth10(:,1)=(5:10:2000)';
    
    ixs = 0.22; ixe = 0.22;  ixd = 0.02; ixw = (1-ixs-ixe-3*ixd)/4;
    iys = 0.04; iye = 0.645;  iyd = 0.05; iyw = (1-iys-iye-0*iyd)/1;
    
    pos{11}  = [ixs+0*ixw+0*ixd   iys+0*iyw+0*iyd  1*ixw 1*iyw]; 
    pos{12}  = [ixs+1*ixw+1*ixd   iys+0*iyw+0*iyd  1*ixw 1*iyw]; 
    pos{13}  = [ixs+2*ixw+2*ixd   iys+0*iyw+0*iyd  1*ixw 1*iyw]; 
    pos{14}  = [ixs+3*ixw+3*ixd   iys+0*iyw+0*iyd  1*ixw 1*iyw]; 
   
    clear color color0 
    color=cbrewer('seq', 'Blues', 60,'pchip');
    color(:,:)=color(60:-1:1,:);

    
% CT. Profiles
subplot('position',pos{11}) 
        for year=1:42
            line00=plot(squeeze(CT_ano_NA(1:end,5,year)),-depth10(1:51,1));
            set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
            hold on
        end
        line00=plot(squeeze(CT_ano_NA(1:end,5,5)),-depth10(1:51,1));
        set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 

        hold on
        line23=plot(squeeze(CT_ano_NA(1:end,5,43)),-depth10(1:51,1));
        set(line23,'color',[0.9,0.3,0.3],'LineWidth',5,'linestyle','-'); 

        text(-1.195,-420,'1981-2022','fontsize',16,'color',color(5,:))
        text(0.6,-420,'2023','fontsize',16,'color',[.9 .3 .3],'FontWeight','bold')

        set(gca,'Ylim',[-500 0],'ycolor','k') 
        set(gca,'YTick',-500:100:0)
        set(gca,'YTickLabel',{'500','400','300','200','100','0'},'fontsize',16)
        set(gca,'Xlim',[-1.2 1.2]) 
        set(gca,'XTick',-1.2:0.6:1.2)
        set(gca,'XTickLabel',{'-1.2',[],'0',[],'1.2'},'fontsize',16)

        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.4,'GridLineStyle','-')
        title('g. Temp. anom. (May)','fontsize',16,'color','k','FontWeight','bold')
        text(0.65,-470,'[ \circC ]','fontsize',16)
        ylabel(['depth [m]'],'fontsize',16)


    
% CT. Profiles
subplot('position',pos{12}) 
        for year=1:42
            line00=plot(squeeze(CT_ano_NA(1:end,6,year)),-depth10(1:51,1));
            set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
            hold on
        end
        line00=plot(squeeze(CT_ano_NA(1:end,6,5)),-depth10(1:51,1));
        set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 

        hold on
        line23=plot(squeeze(CT_ano_NA(1:end,6,43)),-depth10(1:51,1));
        set(line23,'color',[0.9,0.3,0.3],'LineWidth',5,'linestyle','-'); 

        set(gca,'Ylim',[-500 0],'ycolor','k') 
        set(gca,'YTick',-500:100:0)
        set(gca,'YTickLabel',[],'fontsize',16)
        set(gca,'Xlim',[-1.2 1.2]) 
        set(gca,'XTick',-1.2:0.6:1.2)
        set(gca,'XTickLabel',{'-1.2',[],'0',[],'1.2'},'fontsize',16)

        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.4,'GridLineStyle','-')
        title('h. Temp. anom. (June)','fontsize',16,'color','k','FontWeight','bold')
        text(0.65,-470,'[ \circC ]','fontsize',16)
    

% CT. Profiles
subplot('position',pos{13}) 
        for year=1:42
            line00=plot(squeeze(CT_ano_NA(1:end,7,year)),-depth10(1:51,1));
            set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
            hold on
        end
        line00=plot(squeeze(CT_ano_NA(1:end,7,5)),-depth10(1:51,1));
        set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 

        hold on
        line23=plot(squeeze(CT_ano_NA(1:end,7,43)),-depth10(1:51,1));
        set(line23,'color',[0.9,0.3,0.3],'LineWidth',5,'linestyle','-'); 


        set(gca,'Ylim',[-500 0],'ycolor','k') 
        set(gca,'YTick',-500:100:0)
        set(gca,'YTickLabel',[],'fontsize',16)
        set(gca,'Xlim',[-1.2 1.2]) 
        set(gca,'XTick',-1.2:0.6:1.2)
        set(gca,'XTickLabel',{'-1.2',[],'0',[],'1.2'},'fontsize',16)

        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.4,'GridLineStyle','-')
        title('i. Temp. anom. (July)','fontsize',16,'color','k','FontWeight','bold')
        text(0.65,-470,'[ \circC ]','fontsize',16)
    

    
% CT. Profiles
subplot('position',pos{14}) 
        for year=1:42
            line00=plot(squeeze(CT_ano_NA(1:end,8,year)),-depth10(1:51,1));
            set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
            hold on
        end
        line00=plot(squeeze(CT_ano_NA(1:end,8,5)),-depth10(1:51,1));
        set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 

        hold on
        line23=plot(squeeze(CT_ano_NA(1:end,8,43)),-depth10(1:51,1));
        set(line23,'color',[0.9,0.3,0.3],'LineWidth',5,'linestyle','-'); 

        set(gca,'Ylim',[-500 0],'ycolor','k') 
        set(gca,'YTick',-500:100:0)
        set(gca,'YTickLabel',[],'fontsize',16)
        set(gca,'Xlim',[-1.2 1.2]) 
        set(gca,'XTick',-1.2:0.6:1.2)
        set(gca,'XTickLabel',{'-1.2',[],'0',[],'1.2'},'fontsize',16)

        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.4,'GridLineStyle','-')
        title('j. Temp. anom. (August)','fontsize',16,'color','k','FontWeight','bold')
        text(0.65,-470,'[ \circC ]','fontsize',16)

        
% #########################################################################
% #########################################################################
   
disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')


