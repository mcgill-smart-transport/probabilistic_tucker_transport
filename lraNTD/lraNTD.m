function [ Yn,hist] = lraNTD( Y,opts )
%% Nonnegative Tucker Decomposition based on Low-rank Approximation speedup.
%   Usage: Yn = lraNTD_ANLS( Y,opts );
%  Output: Yn is a ttensor
%    opts.
%         NumOfComp: a vector specifying the dimension of Yn.core;
%         nlssolver: [HALS]|APG|MU the solver to be used for nonnegative least-squares
%         corealg: [APG]|MU|none: algorithm for nonnegative core. 'none'
%            will return a real valued core tensor instead of an nonnegative
%            one.
%         maxiter: [100] max number of iterations
%         maxiniter: [20] max number of iterations for internal loop (for each
%               sub-nls problem)
%         tol: 1e-6 the algorithm terminates if ||A(1)-A(1)_old||<tol
%         trackit: [20] check the results after each 'trackit' iterations
%         tdalgFile: a string specifying the file for unconstrained
%               Tucker decomposition. Otherwise 'tucker_als' with 5
%               iterations will be called. 
%         sparsity: an (N+1)-by-1 nonnegative vector whose entries should be
%                significatnly less than 1. sparsity(N+1) is used for the
%                sparsity of the core tensor.
%   This code depends on the TensorToolbox which is available at:
%           http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
%
%  If you think this algorithm is useful, please cite
%     [1]  G. Zhou, A. Cichocki, Q. Zhao, and S. Xie, "Nonnegative Matrix and Tensor Factorizations : An algorithmic perspective," Signal Processing Magazine, IEEE , vol.31, no.3, pp.54,65, May 2014
%          doi: 10.1109/MSP.2014.2298891
%          URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6784087&isnumber=6784017
%     [2]  G. Zhou, A. Cichocki, Q. Zhao, and S. Xie, Efficient Nonnegative Tucker Decompositions: 
%          Algorithms and Uniqueness, [ArXiv e-prints]: http://arxiv.org/abs/1404.4412
%
%   by Guoxu Zhou
%   http://www.bsp.brain.riken.jp/~zhougx/tensor.html
%   E-mail: zhouguoxu@gmail.com or zhouguoxu@ieee.org
%
%   Last updated: April-15, 2014

defopts=struct('NumOfComp',[],'FacAlg','hals','CoreAlg','apg','MaxIter',500,'MaxInIter',20,...
    'MinIter',100,'Tol',1e-3,'TrackIt',20,'TDAlgFile','FileName','L1_Core',0,...
    'L1_Fac',[],'L2_Fac',[],'NN_Core',true,'LRAmode',true,'initFile','FileNameforInit');
if ~exist('opts','var')
    opts=struct();
end
[R,nlssolver,nncore,maxiter,maxiniter,miniter,tol,trackit,tdalgFile,L1_Core,L1_Fac,...
    L2_Fac,NN_Core,LRAmode,initFile]=scanparam(defopts,opts);


I=size(Y);
N=numel(I);

if isempty(L1_Fac)
    L1_Fac=zeros(1,N);
end
if isempty(L2_Fac)
    L2_Fac=zeros(1,N);
end
if min(L1_Fac)<0||min(L2_Fac)<0
    error('L1_Fac and L2_Fac must be nonnegative.');
end
if min(L1_Core)<0
    error('L1_Core must be nonnegative.');
end
if isempty(R)
    error('The parameter ''NumOfComp'' must be specified.');
end
R=round(R);
R=min(R,I);
if any(R<eps)
    error('''NumOfComp'' must be nonzero integers.');
end

if ~NN_Core&&~strcmpi(nncore,'ALS')
    nncore='ALS';
    warning('ALS is applied for the unconstrained core.');
end

activeModes=~isinf(R);
if numel(find(activeModes==1))==1
    error('Please use NMF algorithms directly.');
end

if ~strcmpi(class(Y),'ttensor')
    if LRAmode
        if ~exist(tdalgFile,'file')
            warning('tdalgFile is not specified. MLSVD will be applied to perform unconstrained Tucker decomposition.');
%             [Y]=tucker_als(Y,R,'maxiters',5,'tol',1e-6,'printitn',0);
            [temp,mlsvd_order]=sort(I./R,'descend');
            Y=tensor(Y);
            U=cell(1,N);
            for n=1:N
                currentMode=mlsvd_order(n);
                if activeModes(currentMode)
                    U{currentMode}=nvecs(Y,currentMode,R(currentMode));
                    Y=ttm(Y,U{currentMode},currentMode,'t');
                else
                    U{currentMode}=speye(I(currentMode),I(currentMode));
                end
            end
            Y=ttensor(tensor(Y),U);
            clear U;
        else
            Y=tensor(Y);
            tdalg=load(tdalgFile);
            tdalg.algopts.NumOfComp=R;
            Y=tdalg.alg(Y,tdalg.algopts);
        end
        Y.core=tensor(Y.core);
    else
        Y=tensor(Y);
    end
else
    disp('The data tensor is a ttensor. LRAmode will be adopted automatically.');
    LRAmode=true;
    Y.core=tensor(Y.core);
end
R=round(R);


%% Initialization
A=cell(1,N);
AtA=cell(1,N);
AtyA=cell(1,N);


if exist(initFile,'file')
    init=load(initFile);
    if isfield(init,'A')
        A=init.A;
    else
        A{n}=rand(I(n),R(n));
    end
    if isfield(init,'core')
        core=tensor(init.core);
    else
        core=tensor(rand(R));
    end
    clear init;    
else
    for n=1:N
        if activeModes(n)
            A{n}=rand(I(n),R(n));
        else
            A{n}=speye(I(n),I(n));
        end
    end
    core=tensor(rand(R));
end

for n=1:N
    if activeModes(n)
        AtA{n}=(A{n})'*A{n};
    else
        A{n}=speye(I(n),I(n));
        AtA{n}=A{n};
    end
end



if LRAmode
    AtyA=cell(1,N);
    for n=1:N
        if activeModes(n)
            AtyA{n}=A{n}'*Y.U{n};
        else
            AtyA{n}=speye(I(n),I(n));
        end
    end
end


Xtilde=[];
hist=[];
    if nargout==2
        hist=inf(1,maxiter);
        hist(1)=fitness(Y,ttensor(core,A));
    end
    
n_pos=find(activeModes,1,'last');
activeModesNames=find(activeModes==1);
for it=1:maxiter
    A10=A{1};
    for n=1:N
        if ~activeModes(n)
            continue;
        end
        
        %% update A{n}
         %% compute C=BtB Q=YB
         nindices=activeModesNames(activeModesNames~=n);

         X=ttm(core,AtA(nindices),nindices);
         core_n_T=double(tenmat(core,n,'t'));
         if LRAmode
             Xtilde=ttm(tensor(Y.core),AtyA(nindices),nindices);
         else
             Xtilde=ttm(Y,A(nindices),nindices,'t');
         end
         
         % compute C
         BtB=double(tenmat(X,n));
         BtB=BtB*core_n_T+L2_Fac(n)*eye(R(n),R(n));
         
         % compute Q
         YB=double(tenmat(Xtilde,n));
         YB=YB*core_n_T;
         if LRAmode
            YB=Y.U{n}*YB;
         end
         

             switch upper(nlssolver)
                 case 'HALS'
                     dBtB=max(diag(BtB),eps);
                     %% hals
                     for init=1:maxiniter
%                          od=randperm(R(n));
                         od=1:R(n);
                         for r=od
%                              A{n}(:,r)=max(A{n}(:,r)*dBtB(r)+YB(:,r)-A{n}*BtB(:,r)-L1_Fac(n),eps)/dBtB(r);
                             idx=[1:r-1 r+1:R(n)];
                             A{n}(:,r)=max((YB(:,r)-A{n}(:,idx)*BtB(idx,r)-L1_Fac(n))/dBtB(r),1e-9);
                         end
                     end
                 case 'APG'
                     rho=1/norm(BtB,'fro');
                     alpha1=1;
                     Lk=A{n};
                     for init=1:maxiniter
                         An0=A{n};
                         grad=-YB+Lk*BtB;
                         A{n}=max(Lk-rho*grad-L1_Fac(n),eps);
                         Andiff=A{n}-An0;
                         alpha0=alpha1;
                         alpha1=(1+sqrt(4*alpha0^2+1))/2;
                         Lk=A{n}+((alpha0-1)/alpha1)*Andiff;
                     end
                 case 'MU'
                     YB=YB-L1_Fac(n);
                     for init=1:maxiniter
                         A{n}=A{n}.*(YB)./max(A{n}*BtB,eps);
                         A{n}=max(A{n},eps);
                     end
                 case 'ALS'
                     A{n}=YB/(BtB+(1e-6)*eye(R(n),R(n)));
                     A{n}=max(A{n}-L1_Fac(n),eps);
                 case 'BPP'
                     A{n}=nnlsm_blockpivot(BtB+(1e-6)*eye(R(n),R(n)),YB',1,A{n}');
                     A{n}=max(A{n}-L1_Fac(n),eps)';
                 otherwise
                     error('Unsupported nls-solver.');
             end

         
         nrm=max(abs(A{n}));
         A{n}=bsxfun(@rdivide,A{n},max(nrm,eps));
         core=ttm(core,diag(nrm),n);
         
         AtA{n}=A{n}'*A{n};
         if LRAmode
             AtyA{n}=A{n}'*Y.U{n};
         end
         
    end % for all n     
    
    if rem(it,trackit)==0&&it>miniter
        if max(max(abs(A10-A{1})))<tol
            break;
        end
    end
    
    n=n_pos;
    %% update core
    if NN_Core
        switch upper(nncore)
            case 'ALS'
                iA=cell(1,N);
                [iA([nindices,n])]=cellfun(@(x) inv(x+(1e-6)*eye(size(x))),AtA([nindices,n]),'uni',false);
                
                if LRAmode
                    core=ttm(Xtilde,iA(nindices),nindices);
                    core=ttm(core,iA{n}*(AtyA{n}),n);        
                else
                    core=ttm(Xtilde,iA(nindices),nindices);
                    core=ttm(core,iA{n}*(A{n}'),n);                    
                end
                core=max(double(core),eps);
                core=tensor(core);
            case 'MU'
                if LRAmode
                    enum=double(ttm(Xtilde,AtyA{n},n));
                else
                    enum=double(ttm(Xtilde,A{n},n,'t'));
                end
                enum=enum-L1_Core;
                for init=1:maxiniter
                    core=core.*(enum./max(double(ttm(core,AtA)),eps));                    
                end
            case 'APG'
                rho=1;
                for k=1:N
                    if activeModes(k)==1
                        rho=rho*norm(AtA{k},'fro');
                    else
                        rho=rho*sqrt(I(k));
                    end
                end
                rho=1/rho;
                alpha1=1;
                Lk=core;

                if LRAmode
                    gradY=-ttm(Xtilde,AtyA{n},n);
                else
                    gradY=-ttm(Xtilde,A{n},n,'t');
                end

                for init=1:maxiniter
                    core0=core;
                    grad=ttm(Lk,AtA);
                    grad=double(gradY+grad);
                    core=double(core)-rho*grad;
                    core=max(core-L1_Core,eps);
                    corediff=core-core0;
                    alpha0=alpha1;
                    alpha1=(1+sqrt(4*alpha0^2+1))/2;
                    Lk=tensor(core+((alpha0-1)/alpha1)*corediff);
                end
                core=tensor(core);
            otherwise
                error('Unsupported algorithm for estimating nonnegative core tensor.');
        end
    else
        [iA([nindices,n])]=cellfun(@(x) inv(x+(1e-6)*eye(size(x))),AtA([nindices,n]),'uni',false);
        if LRAmode
            core=ttm(Xtilde,iA(nindices),nindices);
            core=ttm(core,iA{n}*(AtyA{n}),n);        
        else
            core=ttm(Xtilde,iA,nindices);
            core=ttm(core,iA{n}*(A{n}'),n);                    
        end
    end
    
    if nargout==2
        hist(it+1)=fitness(Y,ttensor(core,A));
    end
    
%     if it>50&&~rem(it,trackit)
%         if norm(A{1}-A10,'fro')<tol
%             break;
%         end
%     end
%     fprintf('>>> fitness=%f\n',fitness(Y,ttensor(tensor(core),A)));
%     pause;
end
hist(isinf(hist))=[];

Yn=ttensor(tensor(core),A);

end


