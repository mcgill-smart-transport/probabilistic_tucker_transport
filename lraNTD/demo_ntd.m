%              ==== VERY IMPORTANT ====
% This code needs the support of Tensor Toolbox developed by Tamara Kolda 
%  which is available at:
%         http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
%

I=[20,20,20];
R=[5,6,7];
N=numel(I);

% Generate data;
A=cell(N,1);
for n=1:N
    A{n}=rand(I(n),R(n));
end
Y=ttensor(tensor(rand(R)),A);
Y=tensor(Y);


opts=struct('NumOfComp',R,'FacAlg','hals','MaxIter',100,'MaxInIter',20,...
    'TDAlgFile','call_tucker_als_opts.mat');
tic;
[Ydec]=lraNTD(Y,opts);
toc;
fprintf('Complete. Fit=%f\n',fitness(Y,Ydec));

%%
clc;

I=[20,20,20];
R=[5,6,7];
N=numel(I);

% Generate data;
A=cell(N,1);
for n=1:N
    A{n}=rand(I(n),R(n));
end
Y=ttensor(tensor(rand(R)),A);
Y=tensor(Y);


opts=struct('NumOfComp',R,'FacAlg','hals','MaxIter',100,'MaxInIter',20,...
    'TDAlgFile','call_tucker_als_opts.mat');
tic;
[Ydec]=lraNTD(Y,opts);
toc;
fprintf('Complete. Fit=%f\n',fitness(Y,Ydec));

%%
t = num2cell(data_new,1);
id = sub2ind(cate,t{:});
tten = zeros(cate);
for i = 1:length(id)
    tten(id(i)) = weight(i);
end
%%
