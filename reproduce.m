load('data_all.mat');
addpath('dunson');
addpath('lraNTD');
data = bsxfun(@minus,data,min(data)-1);
data_t = data(data(:,1)>4*3600 & data(:,1)<24.5*3600,:);
data_t(:,1) = data_t(:,1)-4.5*3600;
data_t(:,1) = floor(data_t(:,1)/3600)+1;
cate = max(data_t);
disp(cate);

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
%%
core = [4,3,6,6];
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
load tucker4366.mat;
clc;

percent = [0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100];
res = zeros(size(percent));
for i = 1:length(percent)
    data_s = data_t(randsample(1:length(data_t),round(percent(i)/100*length(data_t))),:);
    temp = num2cell(data_s,1);
    idx = sub2ind(cate,temp{:});
    temp = int_hist(idx,prod(cate))';
    tempid = find(temp>=1);
    temp = temp(temp>=1);
    idx = cell(1,4);
    [idx{:}] = ind2sub(cate,tempid);
    data_new = cell2mat(idx);
    weight = temp;
    
    YY = tucker_EM_aggregate(data_new,weight,core,500,cate,1e-10,false,Ydec_new);
    YY = double(YY);
    raw = tten_raw/sum(tten_raw(:));
    r = (YY-raw).^2;
    res(i) = sum(r(:));
end

%%
% even using 1000 trip
semilogx(percent,res,'-s')


%%
data_s = data_t(randsample(1:length(data_t),round(0.1/100*length(data_t))),:);
temp = num2cell(data_s,1);
idx = sub2ind(cate,temp{:});
temp = int_hist(idx,prod(cate))';
tempid = find(temp>=1);
temp = temp(temp>=1);
idx = cell(1,4);
[idx{:}] = ind2sub(cate,tempid);
data_new = cell2mat(idx);
weight = temp;

YY = tucker_EM_aggregate(data_new,weight,core,500,cate,1e-10,true,Ydec_new);


%%
i = 3;
j = 4;
imagesc(squeeze(sum(sum(double(core),i),j)))