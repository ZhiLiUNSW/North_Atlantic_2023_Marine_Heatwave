%% Main Figure 4: Mixed-Layer Heat Budget Anomaly in 2023, ACCESS-OM2 0.25deg ERA5-Forced Simulation

%% ########################################################################
%% 0. Decomposing SWR penatration into SWR_clim and MLD_clim
%% 0.0 MLT Budget Using Monthly SWR, Qh, and MLD 
%% 0.1 MLT Budget Decomposition: MLD' Term, Using SWR Climatology
%% 0.2 MLT Budget Decomposition: Qsw' Term, Using MLD Climatology
%% 0.3 MLT Budget Decomposition: Qsw_H' Term, Using Monthly MLD and SWR


% #########################################################################
% #########################################################################
%% 1. Plotting: ML Temperature Budget Anomaly in 2023 and Decomposition
%  Bar Charts and Maps ####################################################
% #########################################################################
% #########################################################################
clc; clear
figure('Color',[1 1 1]);  %create a new figure of white color background

ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.045; iye = 0.055;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
pos{3 }  = [ixs+0*ixw+0*ixd    iys+2*iyw+2.5*iyd-0.005  3*ixw+2*ixd 1.05*iyw+0.005]; 

% #########################################################################
clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:); 
color0(6,:)=(color0(6,:)+color0(5,:))./2;   
% #########################################################################

    
% #########################################################################    
% #########################################################################
% a. Bar Charts
   disp(['MLT Decomposition from ACCESS OM2 025...'])
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)  = 60*60*24*31;
      second_2_month(2,1)  = 60*60*24*28.25;
      second_2_month(3,1)  = 60*60*24*31;
      second_2_month(4,1)  = 60*60*24*30;
      second_2_month(5,1)  = 60*60*24*31;
      second_2_month(6,1)  = 60*60*24*30;
      second_2_month(7,1)  = 60*60*24*31;
      second_2_month(8,1)  = 60*60*24*31;
      second_2_month(9,1)  = 60*60*24*30;
      second_2_month(10,1) = 60*60*24*31;
      second_2_month(11,1) = 60*60*24*30;
      second_2_month(12,1) = 60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1980-2010
      count=1;
      for year=1980:2010
          disp(['Mixed layer budget in year#',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;

          
          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Vertical mixing
          mixz0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          mixz0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          mixz(:,:,:,count)=squeeze(nansum(mixz0,4)).*second_2_month; clear mixz0
          mixz(mixz==0)=NaN;
          
          
          % ###############################################################
          % heat flux decomposition
          % (a) SW term 
          Qshortwave0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'swflx');
          % Qshortwave0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat');
          Qshortwave(:,:,:,count)=squeeze(nansum(Qshortwave0,4)).*second_2_month; clear Qshortwave0
          Qshortwave(Qshortwave==0)=NaN;
          
          % (b) LW term
          Qlongwave(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'lw_heat').*second_2_month;

          % (c) Sensible term
          Qsensible(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sens_heat').*second_2_month;

          % (d) Latent term
          Qlatent(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'evap_heat').*second_2_month;
          
          % (e) Shortwave penetration at MLB
          Qswr_mld(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;
          % ###############################################################
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      mixz_clim=nanmean(mixz,4); clear mixz
      
      Qshortwave_clim=nanmean(Qshortwave,4); clear Qshortwave
      Qlongwave_clim=nanmean(Qlongwave,4);   clear Qlongwave
      Qsensible_clim=nanmean(Qsensible,4);   clear Qsensible
      Qlatent_clim=nanmean(Qlatent,4);       clear Qlatent
      Qswr_mld_clim=nanmean(Qswr_mld,4);     clear Qswr_mld
       
      % Removing the shortwave penetration from Qnet and SWR
      Qnet_clim=Qnet_clim+Qswr_mld_clim;
      Qshortwave_clim=Qshortwave_clim+Qswr_mld_clim;
      clear Qswr_mld_clim
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['MLT Budget in Year #2023'])
      % (a) MLT tendency
      dmlt(:,:,1:11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month(:,:,1:11);
      dmlt(:,:,12)=NaN;
      
      % (b) Entrainment
      entrainment(:,:,1:11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month(:,:,1:11);
      entrainment(:,:,12)=NaN;
      
      % (c) Total advection
      adve0(:,:,1:11,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve(:,:,1:11)=squeeze(nansum(adve0,4)).*second_2_month(:,:,1:11); clear adve0
      adve(:,:,12)=NaN;
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,1:11,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet(:,:,1:11)=squeeze(nansum(Qnet0,4)).*second_2_month(:,:,1:11); clear Qnet0
      Qnet(:,:,12)=NaN;
      Qnet(Qnet==0)=NaN;
      

      % (e) Vertical mixing
      mixz0(:,:,1:11,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      mixz0(:,:,1:11,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      mixz(:,:,1:11)=squeeze(nansum(mixz0,4)).*second_2_month(:,:,1:11); clear mixz0
      mixz(:,:,12)=NaN;
      mixz(mixz==0)=NaN;
      
                
      % ###############################################################
      % heat flux decomposition
      % (a) SW term 
      Qshortwave0(:,:,1:11,1)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'swflx');
      Qshortwave(:,:,1:11)=squeeze(nansum(Qshortwave0,4)).*second_2_month(:,:,1:11); clear Qshortwave0
      Qshortwave(:,:,12)=NaN;
      Qshortwave(Qshortwave==0)=NaN;

      % (b) LW term
      Qlongwave(:,:,1:11)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'lw_heat').*second_2_month(:,:,1:11);
      Qlongwave(:,:,12)=NaN;
      
      % (c) Sensible term
      Qsensible(:,:,1:11)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sens_heat').*second_2_month(:,:,1:11);
      Qsensible(:,:,12)=NaN;
      
      % (d) Latent term
      Qlatent(:,:,1:11)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'evap_heat').*second_2_month(:,:,1:11);
      Qlatent(:,:,12)=NaN;
      
      % (e) Shortwave penetration at MLB
      Qswr_mld(:,:,1:11)=ncread(['full_mixed_layer_heat_budget_year_2023.nc'],'sw_heat').*second_2_month(:,:,1:11);
      Qswr_mld(:,:,12)=NaN;
      % ###############################################################
      % Removing the shortwave penetration from Qnet and SWR
      Qnet=Qnet+Qswr_mld;
      Qshortwave=Qshortwave+Qswr_mld;
      clear Qswr_mld
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt        = dmlt-dmlt_clim; 
      entrainment = entrainment-entrainment_clim;
      adve        = adve-adve_clim; 
      Qnet        = Qnet-Qnet_clim; 
      mixz        = mixz-mixz_clim; 

      Qshortwave = Qshortwave - Qshortwave_clim; 
      Qlongwave  = Qlongwave  - Qlongwave_clim; 
      Qsensible  = Qsensible  - Qsensible_clim; 
      Qlatent    = Qlatent    - Qlatent_clim; 
      % ###################################################################

      
    % #########################################################################
      % Basin Mask from ACCESS OM2 0.25 forced by ERA5
      load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
      [Sxy,~,~]=function_Cgrid_Area_Distance((lon_025)',(lat_025)');
      Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
      for month=1:12
          basin_mask_NA0(:,:,month)=basin_mask_NA;
      end
      clear month
      
      dmlt(isnan(basin_mask_NA0))=NaN; 
      entrainment(isnan(basin_mask_NA0))=NaN;
      adve(isnan(basin_mask_NA0))=NaN; 
      Qnet(isnan(basin_mask_NA0))=NaN; 
      mixz(isnan(basin_mask_NA0))=NaN; 

      Qshortwave(isnan(basin_mask_NA0))=NaN; 
      Qlongwave(isnan(basin_mask_NA0))=NaN; 
      Qsensible(isnan(basin_mask_NA0))=NaN; 
      Qlatent(isnan(basin_mask_NA0))=NaN; 
      
      dmlt=dmlt.*Sxy;
      entrainment=entrainment.*Sxy;
      adve=adve.*Sxy;
      Qnet=Qnet.*Sxy;
      mixz=mixz.*Sxy;
      
      Qshortwave=Qshortwave.*Sxy;
      Qlongwave=Qlongwave.*Sxy;
      Qsensible=Qsensible.*Sxy;
      Qlatent=Qlatent.*Sxy;
          
      % NA 0-60N
      dmlt_2023(1:12,1)=squeeze(nansum(nansum(dmlt(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      entrainment_2023(1:12,1)=squeeze(nansum(nansum(entrainment(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      adve_2023(1:12,1)=squeeze(nansum(nansum(adve(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qnet_2023(1:12,1)=squeeze(nansum(nansum(Qnet(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      mixz_2023(1:12,1)=squeeze(nansum(nansum(mixz(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      
      Qshortwave_2023(1:12,1)=squeeze(nansum(nansum(Qshortwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlongwave_2023(1:12,1)=squeeze(nansum(nansum(Qlongwave(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qsensible_2023(1:12,1)=squeeze(nansum(nansum(Qsensible(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
      Qlatent_2023(1:12,1)=squeeze(nansum(nansum(Qlatent(:,1:302,1:12),2),1))./squeeze(nansum(nansum(Sxy(:,1:302),2),1));
    % #########################################################################
      
   
    
    % #########################################################################
      clear bar*
      bar_dmlt(:,1)=dmlt_2023(1:12,1);       % MLT tendency
      bar_dmlt(:,2)=Qnet_2023(1:12,1);       % Qnet
      bar_dmlt(:,3)=Qshortwave_2023(1:12,1); % Qnet - SWR
      bar_dmlt(:,4)=Qlatent_2023(1:12,1);    % Qnet - Latent
      bar_dmlt(:,5)=Qlongwave_2023(1:12,1);  % Qnet - LWR
      bar_dmlt(:,6)=Qsensible_2023(1:12,1);  % Qnet - Sensible
      bar_dmlt(:,7)=entrainment_2023(1:12,1)+mixz_2023(1:12,1);% Entrainment + Vertical mixing
      bar_dmlt(:,8)=adve_2023+(dmlt_2023-Qnet_2023-entrainment_2023-mixz_2023-adve_2023); % Temperature advection + Residual

    % SWR decomposition
      load('plot_2023_Heat_Budget_ACCESS_OM2_025_V7_0_AllTerms_Monthly_MLD_SWR_Qh.mat',...
           'dMLT_Qswr_2023_ano')
      bar_dmlt(:,3)=dMLT_Qswr_2023_ano(1:12,1); % Qnet - SWR
      clear dMLT_Qswr_2023_ano

    subplot('position',pos{3})
       h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.95); 
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
       set(h0(7),'FaceColor',[0.366, 0.574, 0.188],'EdgeColor',[0.366, 0.574, 0.188])
       hold on
       set(h0(8),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])


       % ##################################################################
       % Add SWR decomposition
        hold on
        bar_dmlt_Clim(1:9,1:24)=NaN;
        % Time-Varying MLD' Term
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V7_Decomposition_1_MLD_prime.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,7)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano

        % Time-Varying Qsw' Term
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V7_Decomposition_2_Qsw_prime.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,8)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano

        % Q_sw_H' Term
        load('plot_2023_Heat_Budget_ACCESS_OM2_025_V7_Decomposition_3_Qsw_H_prime.mat',...
             'dMLT_Qswr_2023_ano')
        bar_dmlt_Clim(1:9,9)=dMLT_Qswr_2023_ano(1:9,1); clear dMLT_Qswr_2023_ano   

        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.66); 
           hold on
           set(h0_Clim(7),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(8),'FaceColor',[0.2, 0.3, 0.2],'EdgeColor',[0.2, 0.3, 0.2])
           hold on
           set(h0_Clim(9),'FaceColor',[0.2, 0.8, 0.6],'EdgeColor',[0.2, 0.8, 0.6])

           hold on
           plot((6.10:0.001:6.23),1.4*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.10:0.001:6.23),1.2*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.10:0.001:6.23),1.0*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)  
           hold on
           plot((6.10:0.001:6.23),1.4*ones(length((6.10:0.001:6.23)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.10:0.001:6.23),1.2*ones(length((6.10:0.001:6.23)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.10:0.001:6.23),1.0*ones(length((6.10:0.001:6.23)')),'-','color',[0.2, 0.8, 0.6],'linewidth',4)  

           text(6.25, 1.4,'$\mathbf{MLD^\prime}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
           text(6.25, 1.2,'$\mathbf{{Q_{sw}}^\prime}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
           text(6.25, 1.0,'$\mathbf{{Q_{sw,H}}^\prime}$', 'Interpreter', 'latex','fontsize',17,'FontName', 'Aptos')
       % ##################################################################



        hold on
        legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7) h0(8)],...
               'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Vertical mixing + entrainment','Advection + other minor terms',...
               'Location','northeast','Orientation','vertical','NumColumns',2)
        hold on
        set(legend,'fontsize',17)
        hold on
        legend('boxoff')

        % Legend will show names for each color
        legend() 
        set(gca,'Ylim',[-0.43 1.60],'ycolor','k') 
        set(gca,'YTick',-0.4:0.4:1.6)
        set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2','1.6'},'fontsize',19)
        set(gca,'Xlim',[4.5 8.5]) 
        set(gca,'XTick',4.5:1:8.5)
        set(gca,'XTickLabel',{[],[],[],[]},'fontsize',19)

        text(4.915,-0.54,'May','fontsize',19,'color','k','FontWeight','normal')
        text(5.930,-0.54,'Jun','fontsize',19,'color','k','FontWeight','normal')
        text(6.945,-0.54,'Jul','fontsize',19,'color','k','FontWeight','normal')
        text(7.920,-0.54,'Aug','fontsize',19,'color','k','FontWeight','normal')


        grid on
        set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.2,'GridLineStyle','-')
        ylabel(['[ \circC/month ]'],'fontsize',19,'color','k','FontWeight','normal')

        text(4.56,1.42,'a. MLT budget anomalies','fontsize',19,'color','k','FontWeight','bold')
% #########################################################################    
% #########################################################################
    
    
    
    
% #########################################################################    
% #########################################################################
% Maps for Panels (b-g)
ixs = 0.190; ixe = 0.190;  ixd = 0.020; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.045; iye = 0.055;  iyd = 0.050; iyw = (1-iys-iye-2*iyd)/3;
ixw./iyw

pos{11}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{12}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw];

pos{21}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+2*ixw+2*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 



   % ######################################################################
   disp(['MLT Budget from ACCESS OM2 025...'])
      % ###################################################################
      lon_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','xt_ocean');
      lon_025=lon_025+360;
      lat_025=ncread('full_mixed_layer_heat_budget_year_2023.nc','yt_ocean');
      % second_2_month=60*60*24*30;
      second_2_month(1,1)  = 60*60*24*31;
      second_2_month(2,1)  = 60*60*24*28.25;
      second_2_month(3,1)  = 60*60*24*31;
      second_2_month(4,1)  = 60*60*24*30;
      second_2_month(5,1)  = 60*60*24*31;
      second_2_month(6,1)  = 60*60*24*30;
      second_2_month(7,1)  = 60*60*24*31;
      second_2_month(8,1)  = 60*60*24*31;
      second_2_month(9,1)  = 60*60*24*30;
      second_2_month(10,1) = 60*60*24*31;
      second_2_month(11,1) = 60*60*24*30;
      second_2_month(12,1) = 60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % Heat budget terms during 1980-2010
      count=1;
      for year=1980:2010
          disp(['   MLT Budget in Year #',num2str(year)])
          % (a) MLT tendency
          dmlt(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'mlt_tendency').*second_2_month;

          % (b) Entrainment
          entrainment(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'entrainment').*second_2_month;

          % (c) Total advection
          adve0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_advection');
          adve(:,:,:,count)=squeeze(nansum(adve0,4)).*second_2_month; clear adve0
          adve(adve==0)=NaN;


          % (d) Surafce heat flux
          Qnet0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_sbc');
          Qnet(:,:,:,count)=squeeze(nansum(Qnet0,4)).*second_2_month; clear Qnet0
          Qnet(Qnet==0)=NaN;


          % (e) Shortwave penetration at MLB
          Qswr(:,:,:,count)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sw_heat').*second_2_month;


          % (f) Total mixing
          resd0(:,:,:,1)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_submeso');
          resd0(:,:,:,2)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_k33');
          resd0(:,:,:,3)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_diffusion_temp');
          resd0(:,:,:,4)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'neutral_gm_temp');
          resd0(:,:,:,5)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'residual');% residual=mixdownslope_temp + temp_sigma_diff  
          resd0(:,:,:,6)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_eta_smooth');
          
          resd0(:,:,:,7)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_vdiffuse_diff_cbt');
          resd0(:,:,:,8)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_nonlocal_KPP');
          
          resd0(:,:,:,9)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'temp_rivermix');
          resd0(:,:,:,10)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'sfc_hflux_pme');
          resd0(:,:,:,11)=ncread(['full_mixed_layer_heat_budget_year_',num2str(year),'.nc'],'frazil_3d');
          
          resd(:,:,:,count)=squeeze(nansum(resd0,4)).*second_2_month; clear resd0
          resd(resd==0)=NaN;
          
          count=count+1;
      end
      
      dmlt_clim=nanmean(dmlt,4); clear dmlt
      entrainment_clim=nanmean(entrainment,4); clear entrainment
      adve_clim=nanmean(adve,4); clear adve
      Qnet_clim=nanmean(Qnet,4); clear Qnet
      Qswr_clim=nanmean(Qswr,4); clear Qswr
      resd_clim=nanmean(resd,4); clear resd
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % ###################################################################
      % heat budget terms in 2023
      disp(['Mixed layer budget in year#2023'])
      % (a) MLT tendency
      dmlt(:,:,1:11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','mlt_tendency').*second_2_month(:,:,1:11);
      dmlt(:,:,12)=NaN;
      
      % (b) Entrainment
      entrainment(:,:,1:11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','entrainment').*second_2_month(:,:,1:11);
      entrainment(:,:,12)=NaN;
      
      % (c) Total advection
      adve0(:,:,1:11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_advection');
      adve(:,:,1:11)=adve0(:,:,1:11).*second_2_month(:,:,1:11); clear adve0
      adve(:,:,12)=NaN;
      adve(adve==0)=NaN;
      
      
      % (d) Surafce heat flux
      Qnet0(:,:,1:11,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_sbc');
      Qnet(:,:,1:11)=squeeze(nansum(Qnet0,4)).*second_2_month(:,:,1:11); clear Qnet0
      Qnet(:,:,12)=NaN;
      Qnet(Qnet==0)=NaN;
      
      
      % (e) Shortwave penetration at MLB
      Qswr(:,:,1:11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sw_heat').*second_2_month(:,:,1:11);
      Qswr(:,:,12)=NaN;

      % (f) Total mixing
      resd0(:,:,:,1)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_submeso');
      resd0(:,:,:,2)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_k33');
      resd0(:,:,:,3)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_diffusion_temp');
      resd0(:,:,:,4)=ncread('full_mixed_layer_heat_budget_year_2023.nc','neutral_gm_temp');
      resd0(:,:,:,5)=ncread('full_mixed_layer_heat_budget_year_2023.nc','residual');% residual=mixdownslope_temp + temp_sigma_diff  
      resd0(:,:,:,6)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_eta_smooth');
      
      resd0(:,:,:,7)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_vdiffuse_diff_cbt');
      resd0(:,:,:,8)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_nonlocal_KPP');
      
      resd0(:,:,:,9)=ncread('full_mixed_layer_heat_budget_year_2023.nc','temp_rivermix');
      resd0(:,:,:,10)=ncread('full_mixed_layer_heat_budget_year_2023.nc','sfc_hflux_pme');
      resd0(:,:,:,11)=ncread('full_mixed_layer_heat_budget_year_2023.nc','frazil_3d');
      
      resd(:,:,1:11)=squeeze(nansum(resd0,4)).*second_2_month(:,:,1:11); clear resd0
      resd(:,:,12)=NaN;
      resd(resd==0)=NaN;
      % ###################################################################
      % ###################################################################
      
      
      % ###################################################################
      % Anomalies in 2023
      dmlt=dmlt-dmlt_clim; 
      entrainment=entrainment-entrainment_clim;
      adve=adve-adve_clim; 
      Qnet=Qnet-Qnet_clim; 
      Qswr=Qswr-Qswr_clim; 
      resd=resd-resd_clim; 
      % ###################################################################

     

% #########################################################################       
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. MLT tendency (June 2023)','fontsize',19,'FontWeight','bold')
        
        
         
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,6)'+Qswr(:,:,6)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. Surface flux term','fontsize',19,'FontWeight','bold')

        

subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,6)'+resd(:,:,6)'+entrainment(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. Advection + mixing + entrainment','fontsize',19,'FontWeight','bold')
        

        
% #########################################################################       
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,dmlt(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. MLT tendency (July 2023)','fontsize',19,'FontWeight','bold')
        
        
         
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,Qnet(:,:,7)'+Qswr(:,:,7)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. Surface flux term','fontsize',19,'FontWeight','bold')

        

subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_025,lat_025,adve(:,:,7)'+resd(:,:,7)'+entrainment(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('g. Advection + mixing + entrainment','fontsize',19,'FontWeight','bold')

        
        hBar1 = colorbar('EastOutside','vertical');
        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+3*ixw+2*ixd+0.016 iys+0.5*iyw+0*iyd 0.012 1.0*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',19,'FontName','Arial','LineWidth',1.2,'TickLength',0.060);
        m_text(381,95, '[ \circC/month ]','fontsize',18,'FontWeight','normal')
% #########################################################################
% #########################################################################
   
disp(['>> !!'])
disp(['  Save to left built-in screen...'])
disp(['>> !!'])


