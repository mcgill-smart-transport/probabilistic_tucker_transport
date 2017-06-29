%%
addpath('dunson\');
data = csvread('tensor_trips11042011.csv');
data(:,1) = data(:,1) - 5;
data = data(data(:,1)>0,:);
data = data+1;

[y1,x1] = hist(data(:,2),1:310);
[y2,x2] = hist(data(:,3),1:310);
y = y1+y2;
id = zeros(310,2);
id(:,1) = 1:310;
id(:,2) = y>10000;
id(id(:,2)==1,3) = 1:sum(id(:,2));
id(id(:,2)==0,3) = 0;
disp(length(id));

data(:,2) = id(data(:,2),3);
data(:,3) = id(data(:,3),3);
data = data(data(:,2)>0 & data(:,3)>0,:);
%%
cate = max(data);
data_new = data(randsample(length(data),100000),:);
run = 2000;
burnin = 1000;
thin = 10;
numFactor = 6;
numEffect = 100;
alpha = numFactor/5;
a_alpha = 0.25;
b_alpha = 0.25;

[res,final_pi,final_alpha,final_phi,outseq] = ...
    Gibbs_Dunson_Xing_Sparse(data_new,cate,run,burnin,thin,numFactor,numEffect,alpha,a_alpha,b_alpha,1);

%%
figure;
col = colormap(jet(numFactor));
for i = 1:100
    for f = 1:numFactor
        if res{i}.lambda(f) > 0.04
            plot(res{i}.U{1}(:,f),'color',col(f,:)); hold on;
            disp(f);
        end
    end
end

%% subplots
figure;
numAtt = 3;
subplot(2,3,1);
for i = 1:100
    plot(res{i}.lambda,'s'); hold on;
end
subplot(2,3,2);
hist(outseq{1},100); xlabel(num2str(mean(outseq{1})));
%

subplot(2,3,3);
plot(outseq{2});
subplot(2,3,4);
sim_size = 1;
crav_ten = zeros(numAtt,numAtt);
for idxi = 1:numAtt-1
    for idxj = idxi+1:numAtt
        disp([idxi,idxj]);
        cvt = zeros(numEffect,1);
        temp = zeros(cate(idxi),cate(idxj));
        for eff = 1:numEffect
            u1 = res{eff}.U{idxi};
            u2 = res{eff}.U{idxj};
            lam = res{eff}.lambda;
            for fac = 1:length(lam)
                temp = temp + lam(fac)*(u1(:,fac)*u2(:,fac)');
            end
        end
        w = temp/sum(temp(:));
        w1 = sum(w,1);
        w2 = sum(w,2);
        prd = w2*w1;
        crav_ten(idxi,idxj) = sqrt(sum(sum((w-prd).^2./prd))/min(cate(idxi)-1,cate(idxj)-1));
    end
end
crav_ten = crav_ten+crav_ten';
imagesc(crav_ten);caxis([0,0.6]);

%%
profile on;
res = Gibbs_Dunson_Xing_Sparse(data,cate,10,1,1,numFactor,1,alpha,a_alpha,b_alpha,1);
profile viewer;

%%
%
% cate = max(data);
% temp = num2cell(data,1);
% ind=sub2ind(cate,temp{:});
% counts = zeros(max(data));
% for i = 1:size(ind,1)
%     counts(ind(i)) = counts(ind(i))+1;
% end
% tcount = tensor(counts);

%%
addpath('nonnegfac-matlab-master');

P = ncp(tcount,5,'method','hals','tol',1e-5,'verbose',1,'max_iter',1000);

%%
t1 = reshape(double(tcount),[],1);
t2 = reshape(double(P),[],1);
plot(t1,t2,'s');
w = sum(counts,1);
w  = reshape(w(1,:,:),310,310);
%%
subplot(2,2,1);
plot(P.lambda,'-s')
subplot(2,2,2);
plot(P.U{1});
subplot(2,2,3);
plot(P.U{2});
subplot(2,2,4);
plot(P.U{3});
subplot(2,2,3);
kk = 1;
ww = bsxfun(@times,repmat(P.U{3}(:,kk),1,310),P.U{2}(:,kk)');
imagesc(log(ww));
%%
for kk = 1:5
    figure;
    ww = bsxfun(@times,repmat(P.U{3}(:,kk),1,310),P.U{2}(:,kk)');
    imagesc(log(ww));
    plot(sum(ww)./sum(w)); hold on;
end
%%
min(data)
max(data)
run = 10000;
burnin = 5000;
numFactor = 10;
alpha = 1;
a_alpha = 0.25;
b_alpha = 0.25;
[f_tensor,~,~,~,out] = Gibbs_sample_Dunson_Xing(data, run, burnin, 1,numFactor,1,alpha,a_alpha,b_alpha,1);

%%
plot(out{2})