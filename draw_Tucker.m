clear all;
load('tucker4366.mat');

load('tten_raw.mat');
tucker = Ydec_new;

g = double(tucker);
t = tten_raw/sum(tten_raw(:));

g = squeeze(sum(sum(g,1),2));
t = squeeze(sum(sum(t,1),2));
loglog(g(:),t(:),'s');
%%
load('tucker5388.mat');
tucker = Ydec_new;
g = double(tucker);
t = tten_raw/sum(tten_raw(:));

g = squeeze(sum(sum(g,1),2));
t = squeeze(sum(sum(t,1),2));
figure;loglog(g(:),t(:),'s');
%%
tucker=Ydec_new;
g = double(tucker);
t = tten_raw/sum(tten_raw(:));
g = squeeze(sum(sum(g,1),2));
t = squeeze(sum(sum(t,1),2));
loglog(g(:),t(:),'s');
xlim(10.^[-8,-1]);
axis equal;


% Temporal
g = double(tucker.core);
order = sum(sum(sum(g,4),3),2);
[~,ix1] = sort(order,'descend');

h = fig('units','centimeters','width',6,'height',6,'font','helvetica','fontsize',9,'border','off');
temp = tucker.U{1}(:,ix1);

col = {'red','blue','green','orange','purple'};
mar = {'s','^','o','d','v'};
h = zeros(size(tucker.core,1),1);
for i = 1:size(tucker.core,1)
    h(i) = plot(5:24,temp(:,i),strcat('-',mar{i}),'color',scicol(col{i}),'linewidth',1,...
        'markerfacecolor',scicol(col{i}),'markeredgecolor',scicol(col{i}),'markersize',2); hold on;
end

t = legendflex(h, {'T1','T2','T3','T4'}, 'ref', gca, 'ncol',1, ...
    'fontsize',8,'anchor', {'ne','ne'},'buffer',[-15,-15],'box','off','xscale',0.5);
%c = findall(t,'Type','text'); set(c,'Interpreter', 'latex');

set(gca,'xtick',6:3:24);
set(gca,'ytick',0:0.1:0.4);
xlabel('Time of day \boldmath$t$','interpreter','latex');
ylabel('\boldmath${P}(t)$','interpreter','latex');
xlim([4,25]);
ylim([0,0.4])
set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);


set(gca,'position',[0.20, 0.14, 0.735, 0.77]);

axes('position',[0,0,1,1]);
text(0,1.0,'(a)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');
axis off;

name = 'paper_latex/fig_temporal.pdf';
export_fig(name,'-pdf','-transparent','-nocrop','-rgb');
close all;
system(['start ',name]);

%%


%%
load('tucker_all_5  3  8  8.mat');
new_lambda = cell(1,4);
for i = 1:4
    new_lambda{i} = tucker.U{i};
end
lambda = new_lambda;
save('lambda.mat','lambda');
load('tucker_all_5  3  8  8.mat');

%%
t = cell2mat(pi_weight');
idx = find(t(end,:)==max(t(end,:)));
tucker = lambda{idx};
g = double(tucker);
t = tten_raw/sum(tten_raw(:));
g = squeeze(sum(sum(g,1),2));
t = squeeze(sum(sum(t,1),2));
loglog(g(:),t(:),'s');
xlim(10.^[-8,-1]);
axis equal;

%%
c = double(lambda{2}.core);
ix1 = [2,4,1,5,3];
c = c(ix1,:,:,:);
x = sum(sum(c,3),4)
bsxfun(@rdivide,x,sum(x,1))
sum(x,1)
sum(x,2)
%%

% Temporal
g = double(tucker.core);
order = sum(sum(sum(g,4),3),2);
[~,ix1] = sort(order,'descend');
ix1 = [2,4,1,5,3];

h = fig('units','centimeters','width',6,'height',6,'font','helvetica','fontsize',9,'border','off');
temp = tucker.U{1}(:,ix1);

col = {'red','blue','green','orange','purple'};
mar = {'s','^','o','d','v'};
h = zeros(size(tucker.core,1),1);
for i = 1:size(tucker.core,1)
    h(i) = plot(5:24,temp(:,i),strcat('-',mar{i}),'color',scicol(col{i}),'linewidth',1,...
        'markerfacecolor',scicol(col{i}),'markeredgecolor',scicol(col{i}),'markersize',2); hold on;
end

t = legendflex(h, {'T1','T2','T3','T4','T5'}, 'ref', gca, 'ncol',1, ...
    'fontsize',8,'anchor', {'ne','ne'},'buffer',[-15,-15],'box','off','xscale',0.5);
%c = findall(t,'Type','text'); set(c,'Interpreter', 'latex');

set(gca,'xtick',6:3:24);
set(gca,'ytick',0:0.1:0.5);
xlabel('Time of day \boldmath$t$','interpreter','latex');
ylabel('\boldmath${P}(t)$','interpreter','latex');
xlim([4,25]);
ylim([0,0.5])
set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);


set(gca,'position',[0.20, 0.16, 0.735, 0.77]);

axes('position',[0,0,1,1]);
text(0,1.0,'(a)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');
axis off;
%%
name = 'fig_temporal_5388.pdf';
export_fig(name,'-pdf','-transparent','-nocrop','-rgb');
close all;
system(['start ',name]);

%% passenger
clc;
g = double(tucker.core);
order = sum(sum(sum(g,4),3),1);
[~,ix2] = sort(order,'descend');

h = fig('units','centimeters','width',6,'height',6,'font','helvetica','fontsize',9,'border','off');
temp = tucker.U{2}(:,ix2);

col = {'red','blue','green','yellow'};
mar = {'s','^','o','d'}; %colormap('parula');
h = bar(temp','stack');
set(h,'edgecolor','none');
 set(h(1),'facecolor',scicol('blue'));
 set(h(2),'facecolor',scicol('green'));
% set(h(3),'facecolor',scicol('yellow'));


set(gca,'xtick',1:3);
set(gca,'xticklabel',{'C1','C2','C3'});
set(gca,'ytick',0:0.2:1);
xlabel('Passenger type mode','interpreter','latex');
ylabel('\boldmath${P}(c)$','interpreter','latex');
ylim([0,1]); xlim([0,4]);
set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);


set(gca,'position',[0.20, 0.16, 0.735, 0.77]);
axes('position',[0,0,1,1]);
text(0,1.0,'(b)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');
axis off;
%%
name = 'fig_comp5.pdf';
export_fig(name,'-pdf','-transparent','-nocrop','-rgb');
close all;
system(['start ',name]);

%% Spatial
w = squeeze(sum(sum(double(tucker.core),1),2));
[~,ix3] = sort(sum(w,2),'descend');
[~,ix4] = sort(sum(w,1),'descend');


h = fig('units','centimeters','width',12,'height',5,'font','helvetica','fontsize',7,'border','off');
colormap('parula');
lam = tucker.u{3}(:,ix3);
group = zeros(51,1);
mm = max(lam,[],2);
for i = 1:6
    ff = lam(:,i)==mm;
    group(ff,1) = ix3(i);
    group(ff,2) = sum(lam(ff,:),2);
end
group(:,3) = 1:51;
group1 = sortrows(group,[1,-2]);
idx1 = group1(:,3);

lam = tucker.u{4}(:,ix4);
group = zeros(51,1);
mm = max(lam,[],2);
for i = 1:6
    ff = lam(:,i)==mm;
    group(ff,1) = ix3(i);
    group(ff,2) = sum(lam(ff,:),2);
end
group(:,2) = sum(lam,2);
group(:,3) = 1:51;
group2 = sortrows(group,[1,-2]);
idx2 = group2(:,3);


subplot(1,2,1);
w = double(tucker);
w = squeeze(sum(sum(w,1),2));
w = w(idx1,idx1);
imagesc(w); hold on;% colorbar; set(gca,'clim',[-20,-3]); hold on;
len = 0;
for i = 1:6
    len = len+ sum(group1(:,1)==i);
    plot([0,52],[len+0.5,len+0.5],'k','linewidth',0.5); 
    plot([len+0.5,len+0.5],[0,52],'k','linewidth',0.5); 
end
set(gca,'ydir','normal')
xlabel('Factorization','interpreter','latex');

subplot(1,2,2);
w = tten;
w = squeeze(sum(sum(w,1),2));
w = w/sum(sum(w)); 
w = w(idx1,idx1);
imagesc(w);  hold on;%colorbar;  set(gca,'clim',[-20,-3]); hold on;
len = 0;
for i = 1:6
    len = len+ sum(group1(:,1)==i);
    plot([0,52],[len+0.5,len+0.5],'k','linewidth',0.5); 
    plot([len+0.5,len+0.5],[0,52],'k','linewidth',0.5); 
end
xlabel('Observed','interpreter','latex');
set(gca,'ydir','normal');


%%
clf;
g = double(tucker.core);
g = g(ix1,ix2,ix3,ix4);
gx = squeeze(sum(sum(g,3),4));
imagesc(gx);
sum(gx)
sum(gx,2)
gx = squeeze(sum(sum(g,2),4));
sum(gx)
sum(gx,2)
squeeze(sum(g(4,:,:,:),4))
squeeze(sum(g(4,:,:,:),3))

squeeze(sum(sum(g,1),2));

x = tucker;
x.core(:,1,:,:) = x.core(:,1,:,:);
x.core(:,2,:,:) = x.core(:,1,:,:);
x = squeeze(sum(sum(double(x),3),4));
plot(x/sum(x)); hold on;

x = data_t(data_t(:,2)==3,:); 
x = hist(x(:,1),1:24);x = x/sum(x);
plot(x,'r');
%%
g = double(tucker.core);
g = g(ix1,ix2,ix3,ix4);
h = fig('units','centimeters','width',12,'height',8,'font','helvetica','fontsize',7,'border','off');

colormap('parula');
gx = squeeze(sum(g,2));

for i = 1:5
    subplot(2,3,i);
    imagesc(squeeze(gx(i,:,:)));
    set(gca,'clim',[0,0.06],'xtick',1:6,'xticklabel',{'D1','D2','D3','D4','D5','D6','D7','D8'},'ytick',1:6,...
        'yticklabel',{'O1','O2','O3','O4','O5','O6','O7','O8'});
    text(3.2,0,strcat('T',num2str(i)));
    set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);
end

hh = colorbar('fontsize',5,'color','w');
h_bar = findobj(gcf,'Tag','Colorbar');
set(h_bar,'Location','north');
initpos = get(h_bar,'Position');

pos2 = [0.41   0.51  0.20    0.02];
set(h_bar, 'Position',pos2,'xtick',0:0.03:0.06);



%% for 8 zones
ix1 = [2,4,1,5,3];
ix3 = [4   3 8 6 1 2 7 5];
ix4 = [4     2     3     1     6      7 8     5];
g = double(tucker.core);
g = g(ix1,ix2,ix3,ix4);
h = fig('units','centimeters','width',12,'height',8,'font','helvetica','fontsize',7,'border','off');

colormap('parula');
gx = squeeze(sum(g,2));

for i = 1:5
    subplot(2,3,i);
    imagesc(squeeze(gx(i,:,:)));
    set(gca,'clim',[0,0.06],'xtick',1:8,'xticklabel',{'D1','D2','D3','D4','D5','D6','D7','D8'},'ytick',1:8,...
        'yticklabel',{'O1','O2','O3','O4','O5','O6','O7','O8'});
    text(3.2,0,strcat('T',num2str(i)));
    set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);
end

hh = colorbar('fontsize',5,'color','k');
h_bar = findobj(gcf,'Tag','Colorbar');
set(h_bar,'Location','north');
initpos = get(h_bar,'Position');

pos2 = [0.68   0.41  0.20    0.02];
set(h_bar, 'Position',pos2,'xtick',0:0.03:0.06);
%%
name = 'fig_spatial_t8.pdf';
export_fig(name,'-pdf','-transparent','-nocrop','-rgb');
close all;
system(['start ',name]);

%%
h = fig('units','centimeters','width',6,'height',4,'font','helvetica','fontsize',9,'border','off');

w = double(tucker);
yy = squeeze(sum(sum(sum(w,1),2),4));
yy = yy(idx1);
colormap('parula');
w1 = tucker.u{3}(idx1,ix3);
id = [];
for i = 1:6
    t = find(w1(:,i)==max(w1,[],2));
    [~,it]=sort(yy(t),'descend'); t = t(it);
    id = [id;t];
end
y1 = sum(w1,2);
w1 = w1(id,:);
imagesc(w1'); set(gca,'clim',[0,0.25],'ytick',1:6,'xtick',10:10:50);

xlabel('Zones','interpreter','latex');
ylabel('\boldmath$\theta$(origin)','interpreter','latex');
colormap('parula');

set(gca,'position',[0.20, 0.14, 0.735, 0.77]);
set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);

axes('position',[0,0,1,1]);
text(0,1.0,'(c)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');
axis off;

name = 'paper_latex/fig_origin.pdf';
export_fig(name,'-pdf','-transparent','-nocrop','-rgb');
close all;
system(['start ',name]);

%%
h = fig('units','centimeters','width',6,'height',4,'font','helvetica','fontsize',9,'border','off');

w = double(tucker);
yy = squeeze(sum(sum(sum(w,1),2),4));
yy = yy(idx1);
colormap('parula');
w1 = tucker.u{4}(idx1,ix3);
id = [];
for i = 1:6
    t = find(w1(:,i)==max(w1,[],2));
    [~,it]=sort(yy(t),'descend'); t = t(it);
    id = [id;t];
end
y1 = sum(w1,2);
w1 = w1(id,:);
imagesc(w1'); set(gca,'clim',[0,0.25],'ytick',1:6,'xtick',10:10:50);

xlabel('Zones','interpreter','latex');
ylabel('\boldmath$\theta$(origin)','interpreter','latex');
colormap('parula');

set(gca,'position',[0.20, 0.14, 0.735, 0.77]);
set(gca,'Box','on','TickDir', 'in','TickLength',[.02 .02],'layer','top','color','none',...
    'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','off','XColor',[0,0,0],'YColor',[0,0,0],'ZColor',[0,0,0],'LineWidth',0.5);

axes('position',[0,0,1,1]);
text(0,1.0,'(c)','unit','normalized','VerticalAlignment','top','fontname','times','fontsize',10,'fontweight','bold');
axis off;

name = 'paper_latex/fig_origin.pdf';
export_fig(name,'-pdf','-transparent','-nocrop','-rgb');
close all;
system(['start ',name]);
%%
w = w(ix3,ix4);
imagesc(w); colormap('parula');
ix3 
ix4
subplot(1,2,1);
imagesc(w);
col = cell(6,1);
for i = 1:6
    col{i} = group(group(:,1)==i,3);
end
w = tten;
w = squeeze(sum(sum(w,1),2));
w = w/sum(sum(w)); 
w2 = zeros(6,6);
for i = 1:6
    for j = 1:6
        w2(i,j) = sum(sum(w(col{i},col{j})));
    end
end
subplot(1,2,2);
imagesc(w2);


%%
bar(yy(id)); 
xlim([0.5,51.5]); set(gca,'xtick',5:5:50);
%%
yy = squeeze(sum(sum(sum(w,1),2),3));
yy = yy(idx1);
subplot(2,2,2);
w = double(tucker);
cc = double(tucker.core);
w1 = tucker.u{4}(idx1,ix4);
w1 = w1(id,:);
w1 = bsxfun(@rdivide,w1,sum(w1,2));
imagesc(w1'); set(gca,'clim',[0,1]);
subplot(2,2,4);
bar(yy(id));
xlim([0,51.5])

idx1(id(30))
%%
subplot(1,3,1);
w = double(tucker);
w1 = squeeze(sum(sum(w,1),2));
imagesc(log(w1)); colorbar; set(gca,'clim',[-20,-3]);
set(gca,'ydir','normal')

subplot(1,3,2);
w = tten;
w = squeeze(sum(sum(w,1),2));
w2 = w/sum(sum(w)); 
imagesc(log(w2));  colorbar;  set(gca,'clim',[-20,-3]);
set(gca,'ydir','normal');

subplot(1,3,3);
loglog(reshape(w1,[],1),reshape(w2,[],1),'s');
        

        
%%
g = double(Ydec_new.core);
g = g(ix1,ix2,ix3,ix4);
imagesc(squeeze(g(4,3,:,:)));


%%

