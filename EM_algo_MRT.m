%%
clear all; clc;
dall = load('train_trips11042011.mat');
data = double(dall.trips);
data(:,1) = data(:,1)+ 0.5*data(:,4);
data(:,4) = [];
data = data(data(:,1)>=5*3600 & data(:,1)<=24*3600,:);

[y1,x1] = hist(data(:,2),1:200);
[y2,x2] = hist(data(:,3),1:200);
y = y1+y2;
id = zeros(200,2);
id(:,1) = 1:200
id(:,2) = y>=10000;
id(id(:,2)==1,3) = 1:sum(id(:,2));
id(id(:,2)==0,3) = 0;

data(:,2) = id(data(:,2),3);
data(:,3) = id(data(:,3),3);
data = data(data(:,2)>0 & data(:,3)>0,:);
data(:,1) = data(:,1)-5*3600;
data(:,1) = round(data(:,1)/1800);
data = bsxfun(@minus,data,min(data))+1;
data_new = data(randsample(1:size(data),1000000),:);
%%
profile on;
lca_EM(data_new,4,5000);
profile viewer;
%%
clc;close all;

facBIC = zeros(20,10);
for test = 1:1
    data_new = data(randsample(1:size(data),100000),:);
    parfor fac = 1:20
        [~,~,bic] = lca_EM(data_new,fac,5000,1e-1,false);
        facBIC(fac,test) = bic;
    end
end

%%
plot(facBIC,'-s')

%%
close all;
figure;
fac = 8;
[lambda,pi_weight,loglik,bic] = lca_EM_MAP(data_new,fac,5000,1e-10,true);