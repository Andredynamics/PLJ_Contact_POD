clc; clear; close all;
namefolder = 'Snapshot';
data = loadBasilisk2(namefolder);

%% Mean flow. 

%% U
iii = 100;
    figure('units','centimeters','Position',[5 5 6 12]);
     set(0,'DefaultAxesFontName', 'Times New Roman');
     set(0,'DefaultAxesFontSize',11); 
     
        hold on;box on
        u_mat=real(reshape(data.Umean,data.ny,data.nx)');
        
        u_mat=u_mat/max(abs(u_mat(:)));
        
        contourf(data.y_mat'/data.H,data.x_mat'/data.L,u_mat,...
        linspace(-1,1,100) ...
        ,'LineColor','none'); 
    
    plot(mean(yy1mat,1)/data.H,xx'/data.L,'--w', 'LineWidth', 1.5)
    plot(mean(yy2mat,1)/data.H,xx'/data.L,'--w', 'LineWidth', 1.5) 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);     
    set(gca,'ydir','reverse')
    
   
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');

    
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    set(gca,'xlim',[-.8 0.8])
    title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    
   
     colormap(redblueTecplot(256));

     
     caxis([-1 1]); 

     
     cc = colorbar('location','northoutside');
     cc.Limits=[0 0.3];
     cc.TickLabelInterpreter='latex';
     cc.FontSize=10;
     
 
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

nomefile = 'U_Mean';
cd figure
print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
cd ..
 
 u_mat=real(reshape(data.Um(:,iii),data.ny,data.nx)');

 u_mat=u_mat/max(abs(u_mat(:)));

   
    figure('units','centimeters','Position',[5 5 6 12]);
    hold on;box on
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,u_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);

    set(gca,'ydir','reverse')
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    colormap(redblueTecplot(256));
    set(gca,'xlim',[-.8 0.8])
    caxis([-1 1]);
    cc=colorbar('location','northoutside');
    cc.Limits=[0 0.3];
    cc.TickLabelInterpreter='latex';
    cc.FontSize=10;
     
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile='U_100_Snapshot';
    cd figure
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..



%% V

 v_mat=real(reshape(data.Vmean,data.ny,data.nx))';
 
 v_mat=v_mat/max(abs(v_mat(7:end,:)),[],'all');
 

 figure('units','centimeters','Position',[5 5 6 12]);
        hold on;box on

 
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,v_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 

  
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fmean,data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);
    
  
    set(gca,'ydir','reverse')
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    colormap(redblueTecplot(256));
    set(gca,'xlim',[-.8 0.8])
    caxis([-1 1]);
    cc=colorbar('location','northoutside');
    cc.Limits=[0 0.3];
    cc.TickLabelInterpreter='latex';
    cc.FontSize=10;
        
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile=strcat('V_Mean',namefolder);
    cd figure
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..


    v_mat=real(reshape(data.Vm(:,iii),data.ny,data.nx)');
    
    v_mat=v_mat/max(abs(v_mat(:)));
    
 
    figure('units','centimeters','Position',[5 5 6 12]);
    hold on;box on
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,v_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 

    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);
    set(gca,'ydir','reverse')
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    colormap(redblueTecplot(256));
    set(gca,'xlim',[-.8 0.8])
    caxis([-1 1]);
    cc=colorbar('location','northoutside');
    cc.Limits=[0 0.3];
    cc.TickLabelInterpreter='latex';
    cc.FontSize=10;
     
    
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile='V_100_Snapshot';
    cd figure
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..




%% f

   
    figure('units','centimeters','Position',[5 5 6 12]);
    hold on;box on
    f_mat=real(reshape(data.Fmean,data.ny,data.nx)');
    f_mat=f_mat/max(abs(f_mat(:)));
        
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,f_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none');

    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fmean,data.ny,data.nx)',[0.5 0.5],'--w', 'LineWidth', 1.5);
    
    set(gca,'ydir','reverse')
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    %title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    colormap(redblueTecplot(256));
    set(gca,'xlim',[-.8 0.8])
    caxis([-1 1]);
    cc=colorbar('location','northoutside');
    cc.Limits=[0 0.3];
    cc.TickLabelInterpreter='latex';
    cc.FontSize=10;
     
  
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile='F_Mean_100_Snapshot';
    cd figure 
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..
    
    
    
   
    f_mat=real(reshape(data.Fm(:,iii),data.ny,data.nx)');

    f_mat=f_mat/max(abs(f_mat(:)));
    
    
    figure('units','centimeters','Position',[5 5 8 12]);
    hold on;box on    
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,f_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--w', 'LineWidth', 1.5);
    
    set(gca,'ydir','reverse')
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    set(gca,'XAxisLocation','top');
    %title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    colormap(redblueTecplot(256));
    set(gca,'xlim',[-.8 0.8])
    caxis([-1 1]);
    cc=colorbar('location','eastoutside');
    cc.Limits=[0 0.3];
    cc.TickLabelInterpreter='latex';
    cc.Ticks=[0:0.2:1];
    cc.FontSize=10;
     
     pos = get(gcf,'Position');
     set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
     nomefile='F_100_Snapshot';
     cd figure
     print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
     cd ..
 
close all;clearvars -except xx yy yy1mat yy2mat data namefolder


