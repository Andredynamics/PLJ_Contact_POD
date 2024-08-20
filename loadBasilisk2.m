
%%Func that loads the snapshot from Snapshot folder.

function data = loadBasilisk2(namefolder,tcut)

if nargin==1
    tcut=1.5;
end
cd Snapshot
parameters;

fileID = fopen('time','r');

dataArray = textscan(fileID, '%f%f%[^\n\r]', 'Delimiter', ' ',...
    'MultipleDelimsAsOne', true, 'TextType', 'string',...
    'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);


%data.ttru1=data.ttru1(1:1575);
data.ttru2 = unique(dataArray{:, 2});


DT=0.01;
t = data.ttru1:DT:data.ttru1(end);

%data.ttru1 = [0; data.ttru1];
%data.ttru2 = [0; data.ttru2];


files_dat = dir('t_*.dat');

cd ..


[name,~,~]=natsortfiles({files_dat.name},'\d+\.?\d*'); 

[files_dat.name]=deal(name{:});


n = size(files_dat,1); % number of snapshots


Ustored=[];
Vstored=[];
Fstored=[];
Pstored=[];


 bb = waitbar(0,'Caricamento snapshots.');


for is = 1:n

waitbar(is/n,bb,'Caricamento snapshots.');

cd Snapshot
file1 = importdata(files_dat(is).name);
cd ..

%%file1 = importdata(strcat(namefolder,'/',files_dat(is).name));
strcat(namefolder,'/',files_dat(is).name)


u=file1.data(:,3);             % X velocity component (m/s)
v=file1.data(:,4);             % Y velocity component (m/s)
f=file1.data(:,5);             % Volume fraction 
p=file1.data(:,6);             % Pressure (N/m^2)
%m=sqrt(u.*u+v.*v);             % vel magnitude


Ustored1(is,:)=u';  
Vstored1(is,:)=v'; 
Fstored1(is,:)=f'; 
Pstored1(is,:)=p'; 
%Mstored1(is,:)=m';


% Ustored(is,:)=u';  
% Vstored(is,:)=v'; 
% Fstored(is,:)=f'; 
% Pstored(is,:)=p'; 
% Mstored(is,:)=m';

end


x=file1.data(:,1);             % X coordinate (m)
y=file1.data(:,2);             % Y coordinate (m)


x_mat=reshape(x,ny,nx);
y_mat=reshape(y,ny,nx);


close(bb);

 cc = waitbar(0,'Interpolazione nel tempo.');


for iss = 1:size(Ustored1,2)


  waitbar(iss/size(Ustored1,2),cc,'Interpolazione nel tempo.');

  Ustored(:,iss)= spline(data.ttru1, Ustored1(:,iss),t);
  Vstored(:,iss)= spline(data.ttru1, Vstored1(:,iss),t);
  Fstored(:,iss)= spline(data.ttru1, Fstored1(:,iss),t); 
  Pstored(:,iss)= spline(data.ttru1, Pstored1(:,iss),t);


end

close(cc);



% cut = t>=tcut*L/U;

%cut=t>10*L/U;
%cut=t<9.5*L/U;


% Ustored = Ustored(cut,:);
% Vstored = Vstored(cut,:);
% Fstored = Fstored(cut,:);
% Pstored = Pstored(cut,:);
%Mstored=Mstored(cut,:);




% t = t(cut);
n = length(t);

% data.ttru1=data.ttru1(cut);
% data.ttru2=data.ttru2(cut);




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



ustored = Ustored-U_mean_vec;
vstored = Vstored-V_mean_vec;
fstored = Fstored-F_mean_vec;
pstored = Pstored-P_mean_vec;


ustored=ustored';
vstored=vstored';
fstored=fstored';
pstored=pstored';


U_mean_vec = U_mean_vec';
V_mean_vec = V_mean_vec';
F_mean_vec = F_mean_vec';
P_mean_vec = P_mean_vec';



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
