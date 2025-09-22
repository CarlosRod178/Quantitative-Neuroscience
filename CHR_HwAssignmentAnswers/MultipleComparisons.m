%% Mutliple Comparisons
%% Parts 1 and 2
% The results from the Bonferroni are consistently smaller than those of either 
% the uncorrected form or the Benjamini-Hochberg form, in agreement with its conservative 
% nature. This difference seems to extend as the distance between mu increases.

% Code taken from tutorial values
MU1 = 1;
MU2 = 2;
SIGMA = 2;
alpha = 0.05;

% Get random samples, same n
sampSize = 100;
n = 1000;
pU = ones(1,length(n));
X1 = normrnd(MU1, SIGMA, sampSize, n);
X2 = normrnd(MU2, SIGMA, sampSize, n);

% Calculate values from ttest
[~,pU] = ttest2(X1, X2);

% Determine the proportion of results that fall below desired alpha value
statSigpU = sum(pU<alpha) / length(pU) * 100;
fprintf('percent significant %.2f%% (baseline)', statSigpU)

% Apply the Bonferroni correction
alphaB = alpha/n;
statSigpU_B = sum(pU<alphaB) / length(pU) * 100;
fprintf('percent significant %.2f%% (Bonferroni)', statSigpU_B)

% Apply the Benjamini-Hochberg correction
qBH = alpha;
sortedpU = sort(pU);
critVals = (1:length(sortedpU)) ./n .* qBH;

% Locate the index for the last critical value that is less than or equal
% to its associated ranked p value. Throw sanity check to make sure within
% brackets, then make the new criterion value.
newCritIndex = find(sortedpU <= critVals, 1, 'last');
if isempty(newCritIndex) | newCritIndex <= 0
    % Flags instances when there is no satisfied ranked p less than its
    % associated critical value. 
    disp('No significant values found')
else
    alphaBH = sortedpU(newCritIndex);
    statSigpU_BH = sum(pU<alphaBH) / length(pU) * 100;
    fprintf('percent significant %.2f%% (Benjamini-Hochberg)', statSigpU_BH)
end