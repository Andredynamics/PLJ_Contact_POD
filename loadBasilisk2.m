%%Funzione che carica le snapshot in una matrice. Gli devo passare la
%%cartella e tcut cos'è.

function data = loadBasilisk2(namefolder,tcut)

if nargin==1
    tcut=1.5;
end

%%Mettiamo il comando in una stringa e lo eseguiamo con eval.

cd Snapshot

%%Cosa fa parameters?
parameters;

%%fileID è il puntatore a time per leggere, dove time è la matrice?
%%Chi è time?
fileID = fopen('time','r');

%%Si prende i dati dal file, che si chiama time. Come delimitatore usa lo
%%spazio, credo. Fai prova dopo. 

dataArray = textscan(fileID, '%f%f%[^\n\r]', 'Delimiter', ' ',...
    'MultipleDelimsAsOne', true, 'TextType', 'string',...
    'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);

%Che file è time?
%%Penso contenga gli istanti di tempo, e le iterazioni. 

%%Qua definisce i primi campi dello struct data.

%%Nei primi 2 campi ci mette prima e seconda colonna.  
data.ttru1 = unique(dataArray{:, 2});
%data.ttru1=data.ttru1(1:1575);
data.ttru2 = unique(dataArray{:, 2});

%%
DT=0.01;
t = data.ttru1:DT:data.ttru1(end);

%%Ci concatena uno 0 all'inizio. 
%data.ttru1 = [0; data.ttru1];
%data.ttru2 = [0; data.ttru2];



%% proper files ordering

%%Ok, quindi mette in files_dat tutti i nomi di file che cominciano con
%%NomeCartella\t_*.dat; L'asterisco indica che lì può esserci qualsiasi
%%cosa. 

files_dat = dir('t_*.dat');

cd ..

%%natsortfiles ordina i file in files_dat in ordine alfa-numerico e li schiatta in name. 
[name,~,~]=natsortfiles({files_dat.name},'\d+\.?\d*'); 

%%Li ripassiamo ordinati in files_dat.name; deal è per passare da array
%%(name) a struct (files_dat.name).
[files_dat.name]=deal(name{:});

%%Giustamente, il numero di snapshot è la dimensioni di files_dat. Non
%%dovrebbe esserci .name?
n = size(files_dat,1); % number of snapshots

%%Crea gli array dove andranno U, V, F, P. 
Ustored=[];
Vstored=[];
Fstored=[];
Pstored=[];

%Mstored=[];

%%Fa uscire la barra per aspettare. 
 bb = waitbar(0,'Caricamento snapshots.');

%%Per ogni snapshot:
for is = 1:n
%%Fa passare la barra di prima ad is/n. 
waitbar(is/n,bb,'Caricamento snapshots.');

%%Prendiamo i dati da ogni snapshot, il cui nome è quello contenuto in
%%files_dat.name.

cd Snapshot
file1 = importdata(files_dat(is).name);
cd ..

%%file1 = importdata(strcat(namefolder,'/',files_dat(is).name));
strcat(namefolder,'/',files_dat(is).name)

%%u è la terza colonna, della snapshot, 
%%v è la quarta colonna,
%%f è la quinta, 
%%p è la sesta.

u=file1.data(:,3);             % X velocity component (m/s)
v=file1.data(:,4);             % Y velocity component (m/s)
f=file1.data(:,5);             % Volume fraction 
p=file1.data(:,6);             % Pressure (N/m^2)
%m=sqrt(u.*u+v.*v);             % vel magnitude

%%Passiamo la snapshot in stored, che avrà su ogni riga i valori per un
%%certo istante di tempo. 

Ustored1(is,:)=u';  %% si dovrebbe inizializare ma potrebbero non servire tutti gli snap
Vstored1(is,:)=v'; 
Fstored1(is,:)=f'; 
Pstored1(is,:)=p'; 
%Mstored1(is,:)=m';


% Ustored(is,:)=u';  %% si dovrebbe inizializare ma potrebbero non servire tutti gli snap
% Vstored(is,:)=v'; 
% Fstored(is,:)=f'; 
% Pstored(is,:)=p'; 
% Mstored(is,:)=m';

end

%%Le prime due colonne della snapshot sono le coordinate x ed y. 
x=file1.data(:,1);             % X coordinate (m)
y=file1.data(:,2);             % Y coordinate (m)

%%Rendiamo x ed y matrici. 
x_mat=reshape(x,ny,nx);
y_mat=reshape(y,ny,nx);

%%Chiudiamo la waitbar.
close(bb);

%%Apriamo un'altra waitbar.
 cc = waitbar(0,'Interpolazione nel tempo.');

%%Iteriamo su tutte le snapshot 2 volte?
for iss = 1:size(Ustored1,2)

%%Continuiamo la waitbar.
  waitbar(iss/size(Ustored1,2),cc,'Interpolazione nel tempo.');

%%Andiamo a mettere l'interpolazione di Ustored1 e data.ttru1.
  Ustored(:,iss)= spline(data.ttru1, Ustored1(:,iss),t);
  Vstored(:,iss)= spline(data.ttru1, Vstored1(:,iss),t);
  Fstored(:,iss)= spline(data.ttru1, Fstored1(:,iss),t); 
  Pstored(:,iss)= spline(data.ttru1, Pstored1(:,iss),t);

%%Ma t dove l'ha definito?
end
%%Chiudiamo la waitbar. 
close(cc);


%%Ora se t >= di tcut*L/U, cut vale 1. Chi è tcut?
% cut = t>=tcut*L/U;

%cut=t>10*L/U;
%cut=t<9.5*L/U;

%%Ustored è uguale alla prima riga o alla riga 0 di Ustored, tutte le colonne. Che significa? 
%%No davvero ma che sta facendo, si prende la prima riga soltanto, perchè?
% Ustored = Ustored(cut,:);
% Vstored = Vstored(cut,:);
% Fstored = Fstored(cut,:);
% Pstored = Pstored(cut,:);
%Mstored=Mstored(cut,:);

%%Ma stiamo rovinando tutto co sto cut.

%%Che funzione è t?
% t = t(cut);
n = length(t);

% %%Facciamo la stessa cosa con data.ttru1 e 2;
% data.ttru1=data.ttru1(cut);
% data.ttru2=data.ttru2(cut);



%%Ci troviamo i valori medi delle variabili. 
U_mean_vec = mean(Ustored,1);
V_mean_vec = mean(Vstored,1);
F_mean_vec = mean(Fstored,1);
P_mean_vec = mean(Pstored,1);


%M_mean_vec = mean(Mstored,1);
% 
% U_mean_vec = Ustored(end,:);
% V_mean_vec = Vstored(end,:);
% F_mean_vec = Fstored(end,:);
% P_mean_vec = Pstored(end,:);
% M_mean_vec = Mstored(end,:);


%%Troviamo le variazioni sottraendo la media. 
ustored = Ustored-U_mean_vec;
vstored = Vstored-V_mean_vec;
fstored = Fstored-F_mean_vec;
pstored = Pstored-P_mean_vec;

%%Lo giriamo. 
ustored=ustored';
vstored=vstored';
fstored=fstored';
pstored=pstored';

%%Giriamo anche la media. 
U_mean_vec = U_mean_vec';
V_mean_vec = V_mean_vec';
F_mean_vec = F_mean_vec';
P_mean_vec = P_mean_vec';


%%Schiattiamo tutti i valori nella struct data. 
data.Um=Ustored';
data.Vm=Vstored';
data.Fm=Fstored';
data.Pm=Pstored';

data.Umean = U_mean_vec;
data.Vmean = V_mean_vec;
data.Fmean = F_mean_vec;
data.Pmean = P_mean_vec;

data.um = ustored;
data.vm = vstored;
data.fm = fstored;
data.pm = pstored;

data.rho_l=rho_l;
data.mu_l=mu_l;
data.rho_g=rho_g;
data.mu_g=mu_g;
data.g=g;
data.sigm=sigm;
data.U=U;
data.H=H;
data.L=L ;

data.r_rho=r_rho;
data.r_mu=r_mu;
data.RE=RE;
data.WE=WE;
data.FR=FR;
data.eps=eps;
data.p_rif=p_rif;

data.nx=nx;
data.ny=ny;

%data.Lx_g=Lx_g;
%data.ts=ts;
%data.fs=fs;
data.t=t;
data.dt=DT;


data.x_mat=x_mat;
data.y_mat=y_mat;


%data.Mm=sqrt(data.Um.*data.Um+data.Vm.*data.Vm);
data.n=n;



%data.Rho=data.rho_g+(data.rho_l-data.rho_g)*data.Fm;

end