function [f_list_tensor,final_pi,final_alpha,final_phi,outseq] = ...
    Gibbs_Dunson_Xing_Sparse(data,cate,run,burnin,thin,numFactor,numEffect,alpha,a_alpha,b_alpha,verbose)

% -- Lijun Sun -- %
% -- last modified: April 20, 2015 -- %
% Gibbs_sample_Xing.m implements the gibbs sampling procedures as a function
%
% Please mex int_hist.c and randgamma.c for faster computation
% otherwise use histc, gamrnd from matlab

% Gibbs_sample_Xing.m is adapted from orginal code by
% Zhou, Jing, et al. "Bayesian factorizations of big sparse tensors." arXiv preprint arXiv:1306.1598 (2013),
% which is available at
% https://www.researchgate.net/publication/237077205_Bayesian_factorizations_of_big_sparse_tensors

% f_list_tensor is a list of final estimation for the joint probability distribution
% final_pi/final_alpha/final_phi record the effective samples
% outseq is the sequence of sampling values

% run: total run of MCMC
% burnin: burin of MCMC
% thin: thin of MCMC
% numFactor: number of factors (maximum)
% numEffect: final number of effective samples

% nobs: number of observations
% numAtt: size of attributes
[nobs, numAtt] = size(data);

% initial marginal distribution for categories in each attribute
marg = cell(numAtt,1);
for att = 1:numAtt
    temp = int_hist(data(:,att),cate(att))';
    marg{att} = temp./sum(temp);
end

% alpha0
% alpha = 1;
% hyperpara for alpha (gamma distribution)
% a_alpha = 1;
% b_alpha = 1;
% hyperpara for phi (1)
phi = cell(numAtt,1);
for att = 1:numAtt
    % phi{att} = repmat(reshape(marg{att},[],1),1,numFactor);
    phi{att} = ones(cate(att),numFactor)/cate(att); %uniform prior
    % ra = randgamma(ones(cate(att),numFactor));
    % phi{att} = bsxfun(@rdivide,ra,sum(ra));
end

% prior for phi (dirichlet uniform 1)
prior_phi = cell(numAtt,1);
for att = 1:numAtt
    prior_phi{att} = ones(size(marg{att})); % Dunson prior
    % prior_phi{att} = ones(size(marg{att}))/cate(att); % Perks prior
    % prior_phi{att} =  ones(cate(att),numFactor)/cate(att); % Jeffrey prior
end

Vh = betarnd(1,alpha,[numFactor-1,1]);
pi_weight = zeros(numFactor,1);

logmVh = log(1-Vh);
pi_weight(1:numFactor-1) = Vh.*exp(cumsum([0;logmVh(1:end-1)]));
pi_weight(numFactor) = exp(sum(logmVh));

% record pi_sample and alpha sample for all runs
alpha_sample = zeros(run,1);
pi_sample = zeros(run,numFactor);
% record pi, alpha, phi as final samples
% total size
ttsize = (run-burnin)/thin;
final_pi = zeros(ttsize,numFactor);
final_alpha = zeros(ttsize,1);
final_phi = cell(ttsize,1);

ss = 1;
tempz0 = zeros(nobs,numFactor,numAtt);
for iter = 1:run
    % update latent variable z. 
    % This is the most time consuming part of the code. 
    for att = 1:numAtt
        tempz0(:,:,att) = phi{att}(data(:,att),:);
    end
    probz = bsxfun(@times,prod(tempz0,3),pi_weight');
    probz = bsxfun(@rdivide,probz,sum(probz,2)); % bsxfun is fine as well
    % z = mnrnd(1,probz); (slow) the following code runs faster for
    % generating multinomial distribution
    % z is mask matrix with z(i,j) = 1 if obs i belongs to factor group j
    cump = [zeros(nobs,1),cumsum(probz,2)];
    uf = rand(nobs,1);
    z = zeros(nobs,1);
    for f = 1:numFactor
        z(uf > cump(:,f)) = f;
    end
    
    % sample dirichelet process parameters
    % nh is number of obs belong to each factor group
    % update pi_weight
    nh = int_hist(z,numFactor)';
    rsumnh = sum(nh) - cumsum(nh);
    Vh = randbeta(1+nh, alpha+rsumnh);
    logmVh = log(1-Vh);
    pi_weight = Vh.*exp([0;cumsum(logmVh(1:end-1))]);
    loglastweight = sum(logmVh(1:end-1));
    pi_weight(numFactor) = exp(loglastweight);
    
    
    % sample phi
    for att = 1:numAtt
        tt = hist_expand(data(:,att),z,numFactor,cate(att));% an adapted version of histc(), faster
        gam_mat = bsxfun(@plus,tt,prior_phi{att});
        % diri = gamrnd(gam_mat,1); % uncomment if lightspeed is not available
        diri = randgamma(gam_mat); % using lightspeed matlab package
        diri = bsxfun(@rdivide,diri,sum(diri));
        phi{att} = diri;
    end
    
    % sample alpha
    alpha = gamrnd(a_alpha + numFactor - 1, 1/(b_alpha - loglastweight));
    
    % save samples
    pi_sample(iter,:) = pi_weight';
    alpha_sample(iter) = alpha;
    
    if mod(iter,thin) == 0 && iter > burnin
        final_pi(ss,:) = pi_weight;
        final_alpha(ss) = alpha;
        final_phi{ss} = phi;
        ss = ss+1;
    end
    % print iteration if verbose is on
    if verbose
        if mod(iter,500) == 0, disp(iter); end
    end
end

outseq = cell(2,1);
outseq{1} = alpha_sample;
outseq{2} = pi_sample;

% save final esimation of joint pdf as a tensor
% this part uses Matlab tensor toolbox for efficient expanding the k-tensor
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% otherwise, one can save the result in matrix form and expand it manually
idx = 1;
f_list_tensor = cell(numEffect,1);
for ri = size(final_pi,1)-numEffect+1: size(final_pi,1)
    kt = ktensor(final_pi(ri,:)',final_phi{ri});
    f_list_tensor{idx} = kt;
    idx = idx+1;
end
