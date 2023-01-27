% ============================================================== %
%                    Applied Macro Project                       %
% ============================================================== %

% Author: Sam Allen
% Date  : 17/02/2022 (Revised 20/12/2022)

% This code constructs a bivariate VAR using UK data on unemployment 
% and GDP per capita. (See GDP PC.xls and UR.xls for raw data or 
% project_data.xlsx for the combined speadsheet) 

% The long run identification scheme found in Blanchard and Quah (1989)
% is followed to determine the effects of supply and demand shocks on the 
% economy. 

% The baseline model uses the window 1971-2019 with 4 lags, a VAR(4). 


%% ------------------------------------------------------------ %%
%                     ||  House keeping ||

clear;
clc;
close all;
set(0,'defaulttextinterpreter','latex') % Make text match overleaf

%% ----------------------------------------------------------- %%
                  % Select window and lag length
j1 = 1;
lag_selection = 4 ; 

%  while j1 < 61
%  for window_selector = 1:5           % Uncomment to compare lags
   
window_selector = 1;

% Select shock type
for i1 = 1:2
    shocktype = i1;
    
    
%% ------------------------------------------------------------ %%
% Load data (UK data: period 1971:1 - 2021:4)
data  = xlsread('project_data.xlsx');

X(:,1) = 100*diff(log(data(:,1)));  % Log GDP pc and redefine in % growth  
X(:,2) = data(2:end,2);             % Trim Unemployment 

% Specify sample years
if window_selector     == 1
start_year = 1971;
end_year   = 2019;

elseif window_selector == 2
start_year = 1971;
end_year   = 2007;

elseif window_selector == 3
start_year = 1971;
end_year   = 2021;

elseif window_selector == 4
start_year = 2000;
end_year   = 2019;

elseif window_selector == 5
start_year = 1971;
end_year   = 1987;
end

% Trim data so its between the specified years
X = X((start_year-1971)*4+1:(end_year-1971)*4+3,:);

% Detrend data
X = X - ones(length(X),1)*mean(X); 

%% ------------------------------------------------------------ %%
                          % Build Model: 
                          
lags = lag_selection*2;          % Number of lags
dim  = size(X,2);                % Dimension of X 
% ------------------------------------------------------------- %
% Contruct X Matrix and Y Matrix
LM  = lagmatrix(X,0:lags);
X   = LM(lags+1:end,dim+1:end);  % X Matrix
Y   = LM(lags+1:end,1:dim);      % Y Matrix
obs = size(Y,1);                 % Number of observations
% ------------------------------------------------------------- %
% Regression (X'X)^-1 (X'Y):
beta  = X\Y;
eps   = Y - X*beta;
sigma = cov(eps);
F     = companion(beta(1:end,:)); 
% ------------------------------------------------------------- %
% Find Impact Matrix
LR = inv(eye(size(F))-F);
LR = LR(1:dim,1:dim);
D  = chol(LR*sigma*LR','Lower');

Dp = D;
B  = LR\D;
Fp = F;
Bp = B;
B1 = (B*eps')';
% ------------------------------------------------------------- %
% Point estimate IRF for a shock 
shocks     = zeros(dim,1);

if i1==1
    shocks(shocktype)  = 1;
else 
    shocks(shocktype)  = -1;
end

x          = zeros(dim*lags,1);
x(1:dim,1) = Bp*shocks;
% ------------------------------------------------------------- %
% Generate IRF
T    = 100;

for t = 2:T
    x(:,t) = F*x(:,t-1);
end

yp   = cumsum(x(1,:)');
intp = x(2,:)';
% ------------------------------------------------------------- %
% Bootsrap

for i =1:5000
    
    eb      = eps(ceil(obs*rand(obs,1)),:);
    
    Xn(1,:) = X(1,:);
    Yn(1,:) = Xn(1,:)*beta+eb(1,:);
    
    for t = 2:obs
        Xn(t,:) = [Yn(t-1,:) Xn(t-1,1:end-dim)];
        Yn(t,:) = Xn(t,:)*beta+eb(t,:);
    end
    
    betan  = Xn\Yn;
    en     = Yn-Xn*betan;
    sigman = cov(en);
    F      = companion(betan(1:end,:));
    
    LR     = inv(eye(size(F))-F);
    LR     = LR(1:2,1:2);
    D      = chol(LR*sigman*LR','Lower');
    Dp     = D;
    B      = LR\D;
   
    shocks = zeros(dim,1);
    
    if i1==1
    shocks(shocktype)  = 1;
    else 
    shocks(shocktype) = -1;
    end
    
    x          = zeros(dim*lags,1);
    x(1:dim,1) = B*shocks;
    
    for t = 2:T
        x(:,t) = F*x(:,t-1);
    end
    
    y(:,i)   = cumsum(x(1,:)');
    int(:,i) = x(2,:)';
    
end

conf = 0.95;
x    = (1:T)';
    
yh   = quantile(y',conf+(1-conf)/2)';
y1   = quantile(y',(1-conf)/2)';
ym   = median(y,2);
    
inth = quantile(int',conf+(1-conf)/2)';
int1 = quantile(int',(1-conf)/2)';
intm = median(int,2);

% Attempt at supply and demand
if shocktype == 1 
    inth1 = inth;
    int11 = int1;
    intm1 = intm;
    
    yh1   = yh;
    y11   = y1;
    ym1   = ym;
    
    keepvars = {'inth1','int11','intm1','yh1','y11','ym1','holding','j1','lag_selection','B1','window_selector'};
    clearvars('-except',keepvars{:});
else
end

end

figure;
title('figure 1');
subplot(2,2,1)
[hl,hp] = boundedline(x,ym1 ,[ym1-y11 yh1-ym1], 'alpha','b');
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(hl,'linewidth',1.6);
title('Supply Shock: GDP per Capita','Fontsize',14);
set(gcf,'color','w')
ylabel('Percent deviation') 
xtickangle(45);

subplot(2,2,2)
[jl,jp] = boundedline(x,intm1 ,[intm1-int11 inth1-intm1], 'alpha','b');
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(jl,'linewidth',1.6);
title('Supply Shock: Unemployment ','Fontsize',14);
set(gcf,'color','w') 
xtickangle(45);

subplot(2,2,3)
[hl,hp] = boundedline(x,ym ,[ym-y1 yh-ym], 'alpha','b');
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(hl,'linewidth',1.6);
title('Demand Shock: GDP per Capita','Fontsize',14);
set(gcf,'color','w')
ylabel('Percent deviation') 
xlabel('Time (Quarters)') 
xtickangle(45);

subplot(2,2,4)
[jl,jp] = boundedline(x,intm ,[intm-int1 inth-intm], 'alpha','b');
box on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(jl,'linewidth',1.6);
title('Demand Shock: Unemployment','Fontsize',14);
set(gcf,'color','w')
xlabel('Time (Quarters)') 
xtickangle(45);

%% --------------------------------------------------------------- %%
                    % Structural Shocks

figure; 
title('figure 2');
set(gcf,'color','w')
dates = ones(length(eps),1);

for i = 1:length(eps)
    dates(i,1) = 1971.75 + 0.25*i;
end

v2 = [1973.75 -1.75;1973.75 1.75;1974.25 1.75;1974.25 -1.75;1975.5 -1.75;1975.5 1.75;1975.75 1.75;1975.75 -1.75; ...
    1980.25 -1.75; 1980.25 1.75; 1981.5 1.75; 1981.5 -1.75;1990.5 -1.75;1990.5 1.75; 1991.75 1.75; 1991.75 -1.75; ...
    2008.25 -1.75; 2008.25 1.75; 2009.5 1.75; 2009.5 -1.75];
v3 = [1973.75 -0.6;1973.75 0.6;1974.25 0.6;1974.25 -0.6;1975.5 -0.6;1975.5 0.6;1975.75 0.6;1975.75 -0.6; ...
    1980.25 -0.6; 1980.25 0.6; 1981.5 0.6; 1981.5 -0.6;1990.5 -0.6;1990.5 0.6; 1991.75 0.6; 1991.75 -0.6; ...
    2008.25 -0.6; 2008.25 0.6; 2009.5 0.6; 2009.5 -0.6];
f2 = [1 2 3 4;5 6 7 8;9 10 11 12;13 14 15 16;17 18 19 20];
% ------------------------------------------------------------- %
% Supply Shocks
subplot(2,1,1);
supply=line(dates,B1(:,1),'color','blue'); 
hold on 
patch('Faces',f2,'Vertices',v2,'FaceColor','black','facealpha',.1,'Edgecolor','none');
hold off 
title('Supply Side Shocks','Fontsize',14)
xlim([1972 2020]);
ylim([-1.75 1.75]);
hline = refline(0, 0);
hline.Color = 'k';
xtickangle(45); % change angle of x-axis entries

% Demand Shocks
subplot(2,1,2);
demand=line(dates,B1(:,2),'color','blue'); 
hold on 
patch('Faces',f2,'Vertices',v3,'FaceColor','black','facealpha',.1,'Edgecolor','none');
hold off 
title('Demand Side Shocks','Fontsize',14)
xlim([1972 2020]);
ylim([-0.5 0.5]);
hline = refline(0, 0);
hline.Color = 'k';
xtickangle(45); % change angle of x-axis entries

%% ------------------------------------------------------------ %%
               % Compare windows (start loop at top)

% holding(j1,:)    = inth; 
% holding(j1+1,:)  = inth1;
% holding(j1+2,:)  = int1;
% holding(j1+3,:)  = int11;
% holding(j1+4,:)  = intm;
% holding(j1+5,:)  = intm1;
% holding(j1+6,:)  = yh;
% holding(j1+7,:)  = yh1;
% holding(j1+8,:)  = y1;
% holding(j1+9,:)  = y11;
% holding(j1+10,:) = ym;
% holding(j1+11,:) = ym1;
% 
% j1 = j1+ 12; 
% lag_selection = lag_selection+1; 
% keepvars = {'holding','j1','lag_selection','window_selector'};
%     clearvars('-except',keepvars{:});
% end
% end
% x    = (1:100)';
% 
% % 
% figure;
% title('figure 2');
% subplot(2,2,1)
% plot(x,holding(12,:),'Linewidth',1.6)
% hold on
% plot(x,holding(24,:),'Linewidth',1.6)
% plot(x,holding(36,:),'Linewidth',1.6)
% plot(x,holding(48,:),'Linewidth',1.6)
% plot(x,holding(60,:),'Linewidth',1.6)
% ax = gca;
% ax.XGrid = 'off';
% ax.YGrid = 'on';
% ylim([-0.5 1.5]);
% title('Supply: GDP per Capita','Fontsize',14);
% 
% set(gcf,'color','w')
% ylabel('Percent deviation') 
% xtickangle(45);
% 
% 
% subplot(2,2,2)
% plot(x,holding(6,:),'Linewidth',1.6)
% hold on
% plot(x,holding(18,:),'Linewidth',1.6)
% plot(x,holding(30,:),'Linewidth',1.6)
% plot(x,holding(42,:),'Linewidth',1.6)
% plot(x,holding(54,:),'Linewidth',1.6)
% ax = gca;
% ax.XGrid = 'off';
% ax.YGrid = 'on';
% title('Supply: Unemployment','Fontsize',14);
% set(gcf,'color','w') 
% xtickangle(45);
% 
% 
% subplot(2,2,3)
% plot(x,holding(11,:),'Linewidth',1.6)
% hold on
% plot(x,holding(23,:),'Linewidth',1.6)
% plot(x,holding(35,:),'Linewidth',1.6)
% plot(x,holding(47,:),'Linewidth',1.6)
% plot(x,holding(59,:),'Linewidth',1.6)
% ax = gca;
% ax.XGrid = 'off';
% ax.YGrid = 'on';
% title('Demand: GDP per Capita','Fontsize',14);
% set(gcf,'color','w')
% xlabel('Time (Quarters)') 
% ylabel('Percent deviation') 
% xtickangle(45);
% 
% subplot(2,2,4)
% plot(x,holding(5,:),'Linewidth',1.6)
% hold on
% plot(x,holding(17,:),'Linewidth',1.6)
% plot(x,holding(29,:),'Linewidth',1.6)
% plot(x,holding(41,:),'Linewidth',1.6)
% plot(x,holding(53,:),'Linewidth',1.6)
% ax = gca;
% ax.XGrid = 'off';
% ax.YGrid = 'on';
% legend({'2019-1971','2007-1971','2021-1971','2019-2000','1987-1971'},'Location','southeast');
% title('Demand: Unemployment','Fontsize',14);
% set(gcf,'color','w')
% xlabel('Time (Quarters)') 
% xtickangle(45);
% 
