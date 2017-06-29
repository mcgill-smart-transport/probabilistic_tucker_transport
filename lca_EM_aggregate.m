function [lambda,pi_weight,final_loglik,bic] = ...
    lca_EM_aggregate(data,dataweight,nfactor,maxiter,cate,tol,plot_figure,rawten)

if nargin == 3
    tol = 1e-7;
    plot_figure = false;
    cate = max(data);
end

[nobs,ndim] = size(data);
if plot_figure
    figure; colormap('parula');
end


if nargin ~= 8
    lambda = cell(ndim,1);
    for d = 1:ndim
        t = randgamma(ones(cate(d),nfactor));
        lambda{d} = bsxfun(@rdivide,t,sum(t,1)); %dirichlet uniform prior
    end
    pi_weight = randgamma(ones(1,nfactor));
    pi_weight = pi_weight/sum(pi_weight);
else
    lambda = inilambda;
    pi_weight = iniweight;
end


tempz1 = ones(nobs,nfactor);

final_loglik = nan(maxiter,1);
final_loglik(1) = -Inf;

coverg = 100;
iter = 1;
% begin of EM
while coverg > tol && iter < maxiter
    iter = iter+1;
    for d = 1:ndim
        tt = lambda{d}(data(:,d),:);
        tempz1 = tempz1.*tt;
    end
    probz = bsxfun(@times,tempz1,pi_weight);
    tempz1 = ones(nobs,nfactor);
    Eik = bsxfun(@rdivide,probz,sum(probz,2)); % bsxfun is fine as well
    Eik = bsxfun(@times,Eik,dataweight);
    for d = 1:ndim
        tt = sum_eik(Eik,data(:,d),cate(d));
        lambda{d} = bsxfun(@rdivide,tt,sum(tt,1));
    end
    pi_temp = sum(Eik,1);
    pi_temp = pi_temp/sum(pi_temp);
    pi_weight = pi_temp;
    loglik = sum(dataweight.*log(sum(probz,2)));
    coverg = abs(loglik - final_loglik(iter-1));
    final_loglik(iter) = loglik;
    
    % display iteration
    if mod(iter-1,20) == 0 && plot_figure
        subplot(1,3,1);
        plot(bsxfun(@times,lambda{1},pi_weight),'linewidth',2);
        legend(num2str(pi_weight(:)));
        disp([iter,pi_weight]);
        subplot(1,3,2);
        bar(lambda{2,1}','stack'); ylim([0,1]); xlim([0,nfactor+1]);
        subplot(1,3,3);
        qq = -final_loglik;
        plot(qq,'-s','linewidth',2); xlim([max(iter*0.2,1),iter+1]);
        pause(0.1);
        if mod(iter-1,100) == 0
            [~,idx] = sort(pi_weight,'descend');
            pi_weight = pi_weight(idx);
            for d = 1:ndim
                lambda{d} = lambda{d}(:,idx);
            end
        end
    end
end
disp(iter);
bic = -2*final_loglik(iter)+(nfactor*(sum(cate)-length(cate))+nfactor)*log(nobs);
