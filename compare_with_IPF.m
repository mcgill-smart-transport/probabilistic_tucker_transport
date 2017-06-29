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
load tten_raw;
tten_raw = tten_raw/sum(tten_raw(:));

%%
n = size(data_t,1);
data_n = data_t(randsample(n,round(n*0.01)),:);
temp = num2cell(data_n,1);
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

clc;
Ydec_new = tucker_EM_aggregate(data_new,weight,core,3000,cate,1e-10,true,Ydec);

%%
