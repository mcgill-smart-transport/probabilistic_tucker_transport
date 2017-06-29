clear all;
clc;
data = csvread('freq_tensor_2011.csv');
data = sortrows(data,1);
data_t = data(data(:,2)>4*3600 & data(:,2)<24.5*3600,:);
data_t(:,2) = data_t(:,2)-4.5*3600;
data_t(:,2) = floor(data_t(:,2)/3600)+1;
cate = max(data_t);
disp(cate);
data_t(:,1) = data_t(:,1)+1;
data_t(:,3) = [];

%%
user = max(data_t(:,1));
core = [4,6,6];
Eik = ones(user,prod(core));

iter = iter+1;
tempz0 = ones(nobs,prod(core));
for d = 1:ndim
    tempz0 = times(tempz0,lambda{d}(data(:,d),idx{d}));
end
probz = bsxfun(@times,tempz0,pi_weight);

Eik = bsxfun(@rdivide,probz,sum(probz,2));
Eik = bsxfun(@times,Eik,dataweight);
for d = 1:ndim
    tEik = zeros(nobs,core(d));
    for i = 1:core(d)
        tEik(:,i) = Eik*((idx{d}==i)');
    end
    tt = sum_eik(tEik,data(:,d),cate(d));
    lambda{d} = bsxfun(@rdivide,tt,sum(tt,1));
end

pi_temp = sum(Eik,1);
pi_weight = pi_temp/sum(pi_temp);
loglik = sum(dataweight.*log(sum(probz,2)));
coverg = abs(loglik - final_loglik(iter-1));
final_loglik(iter) = loglik;







