%%
clear all; clc;
create = false;
if create
    tt = 80;
    data1 = csvread('tensor_trips11042011.csv');
    data2 = csvread('tensor_trips12042011.csv');
    data3 = csvread('tensor_trips13042011.csv');
    data4 = csvread('tensor_trips14042011.csv');
    data_all = [data1;data2;data3;data4];
    data = data_all+1;
    [y1,x1] = hist(data(:,3),1:tt);
    [y2,x2] = hist(data(:,4),1:tt);
    y = y1+y2;
    id = zeros(tt,2);
    id(:,1) = 1:tt;
    id(:,2) = y>=10000;
    id(id(:,2)==1,3) = 1:sum(id(:,2));
    id(id(:,2)==0,3) = 0;
    disp(length(id));
    data(:,3) = id(data(:,3),3);
    data(:,4) = id(data(:,4),3);
    data = data(data(:,3)>0 & data(:,4)>0,:);
    save('data_all.mat','data','id');
else
    load('data_all.mat');
end
addpath('dunson');
addpath('lraNTD');
data = bsxfun(@minus,data,min(data)-1);
data_t = data(data(:,1)>4*3600 & data(:,1)<24.5*3600,:);
data_t(:,1) = data_t(:,1)-4.5*3600;
data_t(:,1) = floor(data_t(:,1)/3600)+1;
cate = max(data_t);
disp(cate);

% %%
% fac = 15;
% numberofpara = fac*(sum(cate)-length(cate))+fac;
% data_new = data_t(randsample(1:size(data_t,1),round(0.1*length(data_t))),:);
% [lambda,pi_weight,loglik] = lca_EM(data_new,fac,1000,cate,1e-10,true);
% save(strcat('fac',num2str(fac),'.mat'),'id','lambda','loglik','pi_weight');

%%
clc;
temp = num2cell(data_t,1);
idx = sub2ind(cate,temp{:});
temp = int_hist(idx,prod(cate))';
tempid = find(temp>=5);
temp = temp(temp>=5);
idx = cell(1,4);
[idx{:}] = ind2sub(cate,tempid);
data_new = cell2mat(idx);
weight = temp;

t = num2cell(data_new,1);
id = sub2ind(cate,t{:});
tten = zeros(cate);
for i = 1:length(id)
    tten(id(i)) = weight(i);
end
%% BIC score
clc;close all;
facBIC = zeros(25,1);
for fac = 1:25
    disp(fac);
    [lambda,pi_weight,loglik,bic] = lca_EM_aggregate(data_new,weight,fac,5000,cate,1e-8,false);
    facBIC(fac) = bic;
end
plot(facBIC,'-s')

%%
fac = 15;
numberofpara = fac*(sum(cate)-length(cate))+fac;
data_new = data_t(randsample(1:size(data_t,1),round(0.1*length(data_t))),:);
[lambda,pi_weight,loglik] = lca_EM(data_new,fac,1000,cate,1e-10,true);
save(strcat('fac',num2str(fac),'.mat'),'id','lambda','loglik','pi_weight');

%%
core = [5,3,10,10];
disp(prod(core));
numberpara = prod(core)-1 + sum(core.*(cate-1));
data_new = data_t(randsample(1:size(data_t,1),round(0.1*length(data_t))),:);
[lambda,pi_weight,loglik] = tucker_EM(data_new,core,1000,cate,1e-10,true);
save(strcat('tucker',num2str(fac),'.mat'),'id','lambda','loglik','pi_weight');

%%
profile on; lca_EM(data_new,fac,100,cate,1e-10,false);
profile viewer;
%%
profile on;
tic;lca_EM_MAP(data_new,fac,400,cate,1e-10,false);toc;
profile viewer;
%%
temp = data_t(:,[1,3]);
cate = max(temp);
u = num2cell(temp,1);
idx = sub2ind(cate,u{:});
wm = zeros(cate);
for i = 1:length(idx)
    wm(idx(i)) = wm(idx(i))+1;
end
%%

fac = 15;
addpath('poblano_toolbox_1.1');
ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
ncg_opts.RelFuncTol = 1.0e-20;
ncg_opts.MaxIters = 1000;
ncg_opts.MaxFuncEvals = 15000;
ncg_opts.DisplayIters = 50;

ten = zeros(cate);
temp = num2cell(data_tall,1);
idx = sub2ind(cate,temp{:});
for i = 1:length(idx)
    ten(idx(i)) = ten(idx(i)) +1;
end
res = cp_opt(tensor(ten),fac,'alg_options', ncg_opts);
figure; plot(bsxfun(@times,res.U{1},res.lambda'),'linewidth',2);


inilambda = cell(length(cate),1);
w = zeros(4,fac);
for d = 1:length(cate)
    temp = res.U{d}; temp(temp<0) = 0;
    w(d,:) = sum(temp,1);
    inilambda{d} = bsxfun(@rdivide,temp,sum(temp,1));
end

iniweight = (res.lambda') .* prod(w,1);
iniweight = iniweight/sum(iniweight);

data_new = data_t(randsample(1:size(data_t,1),round(0.2*length(data_t))),:);
[lambda,pi_weight,loglik,bic] = lca_EM(data_new,fac,10000,cate,1e-10,true,inilambda,iniweight);

save(strcat('fac',num2str(fac),'.mat'),'id','lambda','loglik','pi_weight');

%%
fac = 15;
lambda = cell(100,1);
pi_weight = cell(100,1);
loglik = cell(100,1);
bic = cell(100,1);
parfor i = 1:100
    data_new = data_t(randsample(1:size(data_t,1),round(0.2*length(data_t))),:);
    [lambda{i},pi_weight{i},loglik{i},bic{i}] = lca_EM(data_new,fac,10000,cate,1e-5,true);
end
save(strcat('fac',num2str(fac),'.mat'),'id','lambda','loglik','pi_weight');

%%
for i = 1:15
    subplot(4,4,i);
    imagesc((lambda{3}(:,i)*lambda{4}(:,i)'));
end
%% reconstruct
fac = 20;
percent = 1:20;
simsize = 1;
len = length(data_t);
res_fac = zeros(simsize,length(percent));
res_inf = zeros(simsize,length(percent));

temp = num2cell(data_t,1);
idx = sub2ind(cate,temp{:});
raw = zeros(prod(cate),1);
for i = 1:length(idx)
    raw(idx(i)) = raw(idx(i))+1;
end
raw = raw/size(data_t,1);

for per = 1:length(percent)
    disp(per);
    for sim = 1:simsize
        dsize = round(len*percent(per)/100);
        data_new = data_t(randsample(1:len,dsize),:);
        [lambda,pi_weight,loglik,bic] = lca_EM(data_new,fac,5000,cate,1e-2,false);
        
        x = ktensor(pi_weight',lambda);
        tens = reshape(double(x),[],1);
        temp = num2cell(data_new,1);
        idx = sub2ind(cate,temp{:});
        infl = zeros(prod(cate),1);
        for i = 1:length(idx)
            infl(idx(i)) = infl(idx(i))+1;
        end
        infl = infl/sum(infl);
        res_inf(sim,per) = sum((infl-raw).^2);
        res_fac(sim,per) = sum((tens-raw).^2);
    end
end
%%
plot(res_inf,'-s'); hold on;
plot(res_fac,'-ro');


%%

data_new = data_t(randsample(1:size(data_t,1),round(0.2*length(data_t))),:);
ten = zeros(cate);
temp = num2cell(data_new,1);
idx = sub2ind(cate,temp{:});
for i = 1:length(idx)
    ten(idx(i)) = ten(idx(i)) +1;
end
addpath('nonnegfac-matlab-master');

figure;
res = ncp(tensor(ten),15,'method','mu','tol',1e-5,'max_iter',1000,'verbose',2);
plot(res.U{1},'linewidth',2);

%%
data_new = data_t(randsample(1:size(data_t,1),round(length(data_t))),:);
ten = zeros(cate);
temp = num2cell(data_new,1);
idx = sub2ind(cate,temp{:});
for i = 1:length(idx)
    ten(idx(i)) = ten(idx(i)) +1;
end

%%
addpath('poblano_toolbox_1.1');

ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
ncg_opts.RelFuncTol = 1.0e-20;
ncg_opts.MaxIters = 10^3;
ncg_opts.DisplayIters = 20;

res = cp_opt(tensor(ten),10,'alg_options', ncg_opts);
plot(bsxfun(@times,res.U{1},res.lambda'),'linewidth',2);
