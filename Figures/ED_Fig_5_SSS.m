%% Extended Data Figure: Monthly Averaged Surface Salinity Anomalies in 2023 from 1981-2010 Mean
%  May, June, July, August, Monthly Mean 1981-2023
%  Monthly Absolute Salinity from IAP TS, 1981-2023


%% #######################################################################
%% Part 1# Maps for 2023 Monthly SSS Anomaly Relative to 1981-2010 Mean
clc;clear

% #########################################################################
% Monthly SSS during 1981-2010
    count=1;
    time_ann(:,1)=(1981:1:2022);
    for year=1981:2010
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

          SSS_mon_19812010(:,:,1:12,count)=squeeze(SA_depth00);
          clear SA_depth00
          count=count+1;
    end
    clear year count


% 2023 Monthly SSS - IAP 
      load(['SA_depth_IAP_C17_2023.mat'],'SA','lon_IAP','lat_IAP')
      lat_IAP=lat_IAP(:,1);
      SA_depth0(:,:,:)=squeeze(nanmean(SA(1,:,:,1:12),1));
      clear SA
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
          
      SA_depth00(1:340,:,:)=SA_depth0(21:360,:,:);
      SA_depth00(341:360,:,:)=SA_depth0(1:20,:,:); clear SA_depth0

      SSS_mon_2023(:,:,1:12,1)=squeeze(SA_depth00);
      clear SA_depth00 
      

% #########################################################################
      SSS_mon_19812010(SSS_mon_19812010==0)=NaN;
      SSS_mon_2023(SSS_mon_2023==0)=NaN;
      
      % Anomaly is relative to the clim-daily mean
      SSS_clim_19812010=squeeze(nanmean(SSS_mon_19812010(:,:,:,1:30,1),4));
      SSS_ano_2023(:,:,:)=SSS_mon_2023-SSS_clim_19812010;
% #########################################################################



% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.270; ixe = 0.270;  ixd = 0.030; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.040; iye = 0.040;  iyd = 0.0450; iyw = (1-iys-iye-2*iyd)/3;

ixw./iyw

%          [left            bottom      width height]
pos{11}  = [ixs          iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{21}  = [ixs+ixw+ixd  iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{31}  = [ixs          iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{41}  = [ixs+ixw+ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 

% #########################################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:6,:)=color(14:-1:9,:);
color0(7:12,:)=color(7:-1:2,:);
% #########################################################################


subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(SSS_ano_2023(:,:,5)',0,0));
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.6 0.6]);
        title('a. SSS anomaly (May 2023)','fontsize',20,'FontWeight','bold')
        
 

subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(SSS_ano_2023(:,:,6)',0,0));
        shading flat
        hold on

        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.6 0.6]);
        title('b. SSS anomaly (June 2023)','fontsize',20,'FontWeight','bold')
        
        
subplot('position',pos{31})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(SSS_ano_2023(:,:,7)',0,0));
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;

        colormap(gca,color0)
        caxis([-0.6 0.6]);
        title('c. SSS anomaly (July 2023)','fontsize',20,'FontWeight','bold')
        
        
        
subplot('position',pos{41})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,smooth2a(SSS_ano_2023(:,:,8)',0,0));
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        
        colormap(gca,color0)
        caxis([-0.6 0.6]);
        title('d. SSS anomaly (August 2023)','fontsize',20,'FontWeight','bold')
        

        hBar1 = colorbar('EastOutside','vertical');

        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.012 iys+1.5*iyw+0*iyd 0.010 1*iyw+1*iyd]);
        set(hBar1, 'ytick',-0.6:0.1:0.6,'yticklabel',{'<-0.6',[],'-0.4',[],'-0.2',[],'0',[],'0.2',[],'0.4',[],'>0.6'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.052);
        ylabel(hBar1, '[ g kg^{-1} ]','rotation',90);
                
        


%% #######################################################################
%% Part 2# Time Series of 2023 Monthly SSS Anomaly Relative to 1981-2010 Mean
% clc;clear

% % #########################################################################
%     % Basin Mask
%     load('basin_mask_Levitus2012.mat','basin_mask','lat_Lev0','lon_Lev0')
%        basin_mask0(:,:)=basin_mask(:,26:155); clear basin_mask
%        basin_mask0(basin_mask0>1)=NaN;
%        basin_mask0(basin_mask0<1)=NaN;
%        for month=1:12
%            basin_dep_mon(:,:,month)=basin_mask0;
%        end
%        clear month basin_mask0 lon_Lev0 lat_Lev0
%        
%        [Sxy,~,~]=function_Cgrid_Area_Distance((21:1:380)',(-64.5:1:64.5)');
%        Sxy(isnan(basin_dep_mon(:,:,1)))=NaN;
% % #########################################################################
% 
% 
% % #########################################################################
% % Monthly SSS during 1981-2010
%     count=1;
%     time_ann(:,1)=(1981:1:2022);
%     for year=1981:2022
%           disp(['loading SSS... year# ',num2str(year)])
%           load(['SA_depth_IAP_C17_',num2str(year),'.mat'],'SA','lon_IAP','lat_IAP')
%           lat_IAP=lat_IAP(26:155,1);
%           SA_depth0(:,:,:)=squeeze(SA(1,:,26:155,1:12)); % 64.5S-64.5N
%           clear SA
% 
%           SA_depth00(1:340,:,:)=SA_depth0(21:360,:,:);
%           SA_depth00(341:360,:,:)=SA_depth0(1:20,:,:); clear SA_depth0
%           SA_depth00(isnan(basin_dep_mon))=NaN;
% 
%           SA_depth00=SA_depth00.*Sxy;
% 
%           % NA 0-60N
%           SSS_mon_19812022(1:12,count)=squeeze(nansum(nansum(SA_depth00(:,66:125,:),2),1))./squeeze(nansum(nansum(Sxy(:,66:125),2),1));
%           % Eastern and Western NA 0-60N, 40W
%           SSS_mon_19812022_ENA(1:12,count)=squeeze(nansum(nansum(SA_depth00(300:end,66:125,:),2),1))./squeeze(nansum(nansum(Sxy(300:end,66:125),2),1));
%           SSS_mon_19812022_WNA(1:12,count)=squeeze(nansum(nansum(SA_depth00(1:300,66:125,:),2),1))./squeeze(nansum(nansum(Sxy(1:300,66:125),2),1));
%           clear SA_depth00 
%           count=count+1;
%     end
%     clear year count
% 
% 
% % 2023 Monthly SSS
%       load(['SA_depth_IAP_C17_2023.mat'],'SA','lon_IAP','lat_IAP')
%       lat_IAP=lat_IAP(26:155,1);
%       SA_depth0(:,:,:)=squeeze(SA(1,:,26:155,1:12)); % 64.5S-64.5N
%       clear SA
%       
%       SA_depth00(1:340,:,:)=SA_depth0(21:360,:,:);
%       SA_depth00(341:360,:,:)=SA_depth0(1:20,:,:); clear SA_depth0
%       SA_depth00(isnan(basin_dep_mon))=NaN;
%       
%       SA_depth00=SA_depth00.*Sxy;
% 
%       % NA 0-60N
%       SSS_mon_2023(1:12,1)=squeeze(nansum(nansum(SA_depth00(:,66:125,:),2),1))./squeeze(nansum(nansum(Sxy(:,66:125),2),1));
%       SSS_mon_2023_ENA(1:12,1)=squeeze(nansum(nansum(SA_depth00(300:end,66:125,:),2),1))./squeeze(nansum(nansum(Sxy(300:end,66:125),2),1));
%       SSS_mon_2023_WNA(1:12,1)=squeeze(nansum(nansum(SA_depth00(1:300,66:125,:),2),1))./squeeze(nansum(nansum(Sxy(1:300,66:125),2),1));
%       clear SA_depth00 
% 
% 
%       SSS_mon_19812022(SSS_mon_19812022==0)=NaN;
%       SSS_mon_19812022_ENA(SSS_mon_19812022_ENA==0)=NaN;
%       SSS_mon_19812022_WNA(SSS_mon_19812022_WNA==0)=NaN;
%       SSS_mon_2023(SSS_mon_2023==0)=NaN;
%       SSS_mon_2023_ENA(SSS_mon_2023_ENA==0)=NaN;
%       SSS_mon_2023_WNA(SSS_mon_2023_WNA==0)=NaN;
% % #########################################################################
% 
% save plot_2023_SSS_IAP_Monthly_1_Spagatti_V1.mat



% #########################################################################
load plot_2023_SSS_IAP_Monthly_1_Spagatti_V1.mat

%   Anomaly is relative to the clim-daily mean
    SSS_mon_19812010_clim=squeeze(nanmean(SSS_mon_19812022(:,1:30),2));
    SSS_mon_19812010_clim_ENA=squeeze(nanmean(SSS_mon_19812022_ENA(:,1:30),2));
    SSS_mon_19812010_clim_WNA=squeeze(nanmean(SSS_mon_19812022_WNA(:,1:30),2));
    count=1;
    for year=1:size(SSS_mon_19812022,2)
        SSS_mon_19812022_ano(:,count)=SSS_mon_19812022(:,count)-SSS_mon_19812010_clim;
        SSS_mon_19812022_ano_ENA(:,count)=SSS_mon_19812022_ENA(:,count)-SSS_mon_19812010_clim_ENA;
        SSS_mon_19812022_ano_WNA(:,count)=SSS_mon_19812022_WNA(:,count)-SSS_mon_19812010_clim_WNA;
        count=count+1;
    end
    
    SSS_mon_2023_ano(:,1)=SSS_mon_2023(:,1)-SSS_mon_19812010_clim;
    SSS_mon_2023_ano_ENA(:,1)=SSS_mon_2023_ENA(:,1)-SSS_mon_19812010_clim_ENA;
    SSS_mon_2023_ano_WNA(:,1)=SSS_mon_2023_WNA(:,1)-SSS_mon_19812010_clim_WNA;    
    
    
% #########################################################################
% #########################################################################
% Plotting ################################################################
% Twitter part
clc;
clear mld* basin*
pos{51}  = [ixs+0.0*ixw  iys+0*iyw+0*iyd   2.0*ixw+1*ixd 0.918*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
subplot('position',pos{51})
    for year=1:42
        line00=plot(1:12,SSS_mon_19812022_ano(:,year));
        set(line00,'color',color(year,:),'LineWidth',2,'linestyle','-'); 
        hold on
    end
    line00=plot(1:12,SSS_mon_19812022_ano(:,5));
    set(line00,'color',color(5,:),'LineWidth',2,'linestyle','-'); 
    
    
% #########################################################################
    % 2023 SSS anomaly from IAP
    hold on
    line23=plot(1:12,SSS_mon_2023_ano(:,1));
    set(line23,'color',[0.950,0.325,0.098],'LineWidth',6,'linestyle','-'); 
% #########################################################################

    
    leg=legend([line00 line23],'1981-2022','2023','Location','northwest');
    set(leg,'fontsize',20)
    hold on
    legend('boxoff')
    
    
    set(gca,'Ylim',[-0.15 0.2],'ycolor','k') 
    set(gca,'YTick',-0.3:0.1:0.3)
    set(gca,'YTickLabel',{'-0.3','-0.2','-0.1','0','0.1','0.2','0.3'},'fontsize',20)
    set(gca,'Xlim',[1 12]) 
    set(gca,'XTick',1:1:12)
    % set(gca,'XTickLabel',{'Jan','Mar','May','Jul','Sep','Nov'},'fontsize',24)
    set(gca,'XTickLabel',{'Jan.','Feb.','Mar.','Apr.','May','Jun.',...
                          'Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'},'fontsize',20)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ g kg^{-1} ]'],'fontsize',20,'color','k','FontWeight','normal')

    title('e. North Atlantic SSS anomaly: 1981-2023','fontsize',20,'color','k','FontWeight','bold')

 
        
% #########################################################################     

disp(' >>')
disp(' Save to the built-in screen')
disp(' >>')



