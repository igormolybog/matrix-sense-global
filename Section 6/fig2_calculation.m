% Initial setup

N = 9; % (N might be equal to 9, 14 or 30)
Rs = 0:0.25:10;

I = size(Rs,2);
results_ball = cell(I,1); % results of optimization of LMI_P will be returned here
results_sphere = cell(I,1); % results of optimization of LMI_P_dB will be returned here
omega = ones(N,1);   % IMPORTANT TO HAVE A COLUMN-VECTOR

%%
% optimization (this might take time, use parfor if so, but be aware
% of memory constraints)

for i = 1:I
    tic
    [results_ball{i}, results_sphere{i}] = BO_func(N, Rs(i), omega);
    toc
end