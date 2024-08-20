clc; clear; close all;

%%Nome della cartella. 
namefolder = 'Snapshot';

%%Schiattiamo le snapshot e i parametri in data.
data = loadBasilisk2(namefolder);

%% Mean flow. 

%%Ora che si fa?

%% U
iii = 100;
%%Apriamo una figure. 
    figure('units','centimeters','Position',[5 5 6 12]);
     set(0,'DefaultAxesFontName', 'Times New Roman');
     set(0,'DefaultAxesFontSize',11); 
     
        hold on;box on

        %%Schiattiamo in u_mat la Umedia, con ny ed nx righe e colonne.
        u_mat=real(reshape(data.Umean,data.ny,data.nx)');
        
        %%Adimensionalizziamo tutto con il massimo.
        u_mat=u_mat/max(abs(u_mat(:)));
        %%Fa il contour riempito di y adimensionalizzato, x adimensionalizzato e la
        %%u_mat.
        contourf(data.y_mat'/data.H,data.x_mat'/data.L,u_mat,...
        linspace(-1,1,100) ...
        ,'LineColor','none'); 
    %%Plotta la media della prima riga yy1mat ed yy2mat, con xx. 
    %plot(mean(yy1mat,1)/data.H,xx'/data.L,'--w', 'LineWidth', 1.5)
    %plot(mean(yy2mat,1)/data.H,xx'/data.L,'--w', 'LineWidth', 1.5) 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);
    %%Delle options.       
    set(gca,'ydir','reverse')
    
    %%Nomi assi. 
    ylabel('$x/L$','interpreter','latex');
    xlabel('$y/H$','interpreter','latex');

    %%Altre options. 
    set(gca,'ticklabelinterpreter','latex','ytick',0:0.2:1,'xtick',-1:0.5:1)
    set(gca,'xlim',[-.8 0.8])
    %title(['$f\approx' num2str(round(fuvC(ind_f)),'%.2f') '$ [Hz]'],'interpreter','latex')
    
     %%Fissiamo la colormap. 
     colormap(redblueTecplot(256));

     %%Fissiamo max e min dei colori. 
     caxis([-1 1]); 

     %%Facciamo la barra dei colori. 
     cc = colorbar('location','northoutside');
     cc.Limits=[0 0.3];
     cc.TickLabelInterpreter='latex';
     cc.FontSize=10;
     
%%Ci da la position di gcf. 
pos = get(gcf,'Position');
%%Altre options. 
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%%Crea nome del file e ci salva l'immagine.  
nomefile = 'U_Mean';
cd figure
print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
cd ..

%% Vediamo il contact alla 200-esima snapshot. 
%%Mettiamo in u_mat Um, tutte le righe, solo 100-esima colonna. 
 u_mat=real(reshape(data.Um(:,iii),data.ny,data.nx)');
 %%Lo adimensionalizziamo col valore max. 
 u_mat=u_mat/max(abs(u_mat(:)));

    %%Creiamo la figure.
    figure('units','centimeters','Position',[5 5 6 12]);
    hold on;box on
    %%Facciamo il contour riempito di questa riga?
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,u_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 
    %%Lo facciamo anche per la 100-esima riga di f. 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);
    
    %%Stesse cose di prima. 
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
    nomefile='U_100_Snapshot';
    cd figure
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..



%% V

 v_mat=real(reshape(data.Vmean,data.ny,data.nx))';
 
 v_mat=v_mat/max(abs(v_mat(7:end,:)),[],'all');
 
 %%Creiamo la figure. 
 figure('units','centimeters','Position',[5 5 6 12]);
        hold on;box on

 %%Facciamo il contour riempito di v. 
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,v_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 

    %%e di f. 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fmean,data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);
    
    %%Alcune options. 
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
      
    %%Salviamo la figura.   
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile=strcat('V_Mean',namefolder);
    cd figure
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..

    %%Lo stiamo rifacendo ma per Vm. 
    v_mat=real(reshape(data.Vm(:,iii),data.ny,data.nx)');
    
    v_mat=v_mat/max(abs(v_mat(:)));
    
    %%Apriamo la figure. 
    figure('units','centimeters','Position',[5 5 6 12]);
    hold on;box on
    %%Facciamo il contour colorato di v. 
    contourf(data.y_mat'/data.H,data.x_mat'/data.L,v_mat,...
    linspace(-1,1,100) ...
    ,'LineColor','none'); 

    %%Contour è per la linea bianca. 
    contour(data.y_mat'/data.H,data.x_mat'/data.L,reshape(data.Fm(:,iii),data.ny,data.nx)',[0.5 0.5],'--b', 'LineWidth', 1.5);
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
     
    %%Salviamo la figure. 
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile='V_100_Snapshot';
    cd figure
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..




%% f
    %%Stessa cosa ma con f. 

    %%Apriamo la figure. 
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
     
    %% Salviamo la figure. 
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    nomefile='F_Mean_100_Snapshot';
    cd figure 
    print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
    cd ..
    
    %%Come prima, lo rifacciamo con f, prima era la media, ma solo per la
    %%100-esima colonna. 
    
    %%Plotta della 100-esima snapshot la f?
    f_mat=real(reshape(data.Fm(:,iii),data.ny,data.nx)');

    f_mat=f_mat/max(abs(f_mat(:)));
    
    %%Apriamo la figure. 
    figure('units','centimeters','Position',[5 5 8 12]);
    hold on;box on    
    %%Fa il contour sia di f_mat che di Fm.  
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
     
     %%Salviamo la figure. 
     pos = get(gcf,'Position');
     set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
     nomefile='F_100_Snapshot';
     cd figure
     print(gcf,nomefile,'-dpdf','-r300'); savefig(strcat(nomefile,'.fig'));
     cd ..


%%Cleariamo le variabili.  
close all;clearvars -except xx yy yy1mat yy2mat data namefolder


