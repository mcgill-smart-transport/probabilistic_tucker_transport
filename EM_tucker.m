%%
clear all;
clc;
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
    id(:,2) = y>=0;
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
addpath('C:\Users\Lijun\Dropbox (¸öÈË)\CommonMatlab\tensor_toolbox_2.6');
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
tempid = find(temp>=1);
temp = temp(temp>=1);
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

% idx = find(weight>=50);
% weight = weight(idx);
% data_new = data_new(idx,:);


%%
cores = [4,3,5,5;
    4,3,6,6;
    4,3,8,8;
    5,3,6,6;
    5,3,8,8];

for kk = 1:5
    core = cores(kk,:);
    disp(prod(core));
    numberpara = prod(core)-1 + sum(core.*(cate-1));
    disp(numberpara);
    
    % APG, MU, HALS
    opts=struct('NumOfComp',core,'FacAlg','HALS','MaxIter',1000,'MaxInIter',25,...
        'TDAlgFile','call_tucker_als_opts.mat');
    
    [Ydec]=lraNTD(tten,opts);
    figure;
    subplot(1,2,1); plot(bsxfun(@rdivide,Ydec.U{1},sum(Ydec.U{1})),'linewidth',3);
    subplot(1,2,2); bar(bsxfun(@rdivide,Ydec.U{2},sum(Ydec.U{2}))','stack');
    
    %%
    % clc;
    % Ydec_new = tucker_EM_aggregate(data_new,weight,core,5000,cate,1e-10,true,Ydec);
    % sum(sum(double(Ydec_new.core),4),3);
    % w = squeeze(sum(sum(double(Ydec_new.core),1),2));
    % w1 = sum(w,1); w2 = sum(w,2);
    % [~,id1] = sort(w1); [~,id2] = sort(w2);
    % w = w(id2,id1);
    % imagesc(log(w));
    % lambda = Ydec_new.U;
    % save('tucker.mat','id','lambda','Ydec_new');
    %plot(reshape(double(Ydec),[],1),reshape(double(Ydec_new),[],1),'s');
    
    % %%
    % fac = 15;
    % numberofpara = fac*(sum(cate)-length(cate))+fac;
    % disp(numberofpara);
    % lca_EM_aggregate(data_new,weight,fac,5000,cate,1e-10,true);
    
    %%
    clc;
    run_size = 10;
    lambda = cell(run_size,1);
    pi_weight = cell(run_size,1);
    loglik = cell(run_size,1);
    bic = cell(run_size,1);
    
    parfor i = 1:run_size
        opts=struct('NumOfComp',core,'FacAlg','HALS','MaxIter',1000,'init','random','MaxInIter',25,...
            'TDAlgFile','call_tucker_als_opts.mat');
        [Ydec]=lraNTD(tten,opts);
        [lambda{i},pi_weight{i},loglik{i}] = tucker_EM_aggregate(data_new,weight,core,2000,cate,1e-5,false,Ydec);
    end
    save(strcat('tucker_all_',num2str(core),'.mat'),'id','lambda','loglik','pi_weight');
    
end

% tucker_EM_aggregate(data_new,weight,core,10000,cate,1e-10,true);
% %%
% profile on;
% [lambda,pi_weight,loglik] = tucker_EM_aggregate(data_new,weight,core,5,cate,1e-10,true);
% profile viewer;