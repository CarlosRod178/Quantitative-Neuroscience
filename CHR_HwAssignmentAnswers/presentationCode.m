%%% Code for QNC Presentation 1
% 
% Start by setting some base variables like number of parcellations used in
% Gordon 2016 atlas. Then load in relevant data for participant
numParcels = 333; 

% Start by loading and transposing the resting state fMRI residual data.
% Then re-sort based on anatomic locations and groupings, using the index2 
% for clarity later on
errtsUnsorted = load('003/003.errts.mni.parcels.1d')';
indSort = load('003/all_sorted.indices2.1d');
errts = zeros(size(errtsUnsorted));
errts(:) = errtsUnsorted(indSort,:);

% Locate all the files in the appropriate directory containing BOLD fMRI
% coefficients. Load them into a holder matrix and then average
coeffsLocate = dir(fullfile('003/', '*_npu.tms.coef.mni.parcels.1d'));
coeffsHolder = ones(length(coeffsLocate), numParcels);
for i = 1:length(coeffsLocate)
    coeffsHolder(i,:) = load(['003/', coeffsLocate(i).name]);
end
coeffsUnsorted = mean(coeffsHolder, 1);
coeffs = zeros(size(coeffsUnsorted));
coeffs(:) = coeffsUnsorted(indSort);
% coeffsUnsorted = coeffsHolder(1,:);
% Load the E-Field data from SIMNIBS and re-sort
eFieldUnsorted = load('003/simnibs_003_srt-mnt.rdlpfc_MRI-B91_A.nii.gz_vector_23/subject_volumes/003_TMS_1-0001_MRI-B91_A_nii_scalar_normE.mni.parcels.1d'); 
eField = zeros(size(eFieldUnsorted));
eField(:) = eFieldUnsorted(indSort);

%% Assess the rsfMRI residuals 
%
% These are the series of BOLD signals across parcels

% Visualize the fMRI residuals as a grayscale imagesc. Note, the
% parcellations are *not* in ascending order since the re-sort
figure()
errtsFig = imagesc(errts(:,:));
colormap('gray');
xlabel('Timeseries')
ylabel('Parcel')
yticklabels(indSort)
title('Resting State fMRI Residuals')

% Plot random residuals against one another to assess if really linear and
% warrants Pearson's R statistic
subNum = 9;
figure()

for i = 1:subNum
    subplot(3,3, i)
    parcelChoices = randi(numParcels, [1,2]);
    xVal = errts(parcelChoices(1),:);
    yVal = errts(parcelChoices(2),:);
    scatPlot = plot(xVal, yVal, '.', MarkerSize=10);
    test = polyfit(xVal, yVal, 1);
    refline(test)
    xlabel(['Parcel ', int2str(parcelChoices(1))]);
    ylabel(['Parcel ', int2str(parcelChoices(2))]);
end

%% Compute Functional Connectivity Matrix

% Calculate the connectivity matrix using a Pearson's correlation matrix.
% This is the standard in the field. All NaN values sent to zero
connMatrix = corr(errts', 'Type', 'Pearson');
connMatrix(isnan(connMatrix)) = 0;

% Display the connectivity using a standard imagesc
figure()
imagesc(connMatrix)
xlabel('Parcel')
ylabel('Parcel')
title('Functional Connectivity Matrix')

% Calculate the BOLD predictions using just the resting state fMRI
predicted = ones(1,333) * connMatrix;
deltaPredicted = eField * connMatrix; 

% Calculate the partial linear correlation when controlling for variables
% in the resting state fMRI. Gives us a sense of how much of these results
% are due to the E-Field data
partialCorr = partialcorr(coeffs', deltaPredicted', predicted');
pairwiseR = corr(coeffs', deltaPredicted', Rows="Pairwise");
pairwiseRNoEF = corr(coeffs', predicted', Rows="Pairwise");

figure()
plot(coeffs, deltaPredicted, '.', MarkerSize=10)
xlabel('TMS BOLD Response')
ylabel('Predicted Response')
title(['Partial Coefficient ', num2str(partialCorr)])

%% ChatGPT Laplacian Distribution for visualization
% Parameters
mu = 0;       % mean (center)
b = 1;        % scale (controls width)

% Range of x values
x = -10:0.01:10;

% Laplacian PDF
f = (1/(2*b)) * exp(-abs(x - mu)/b);

% Plot
figure;
plot(x, f, 'LineWidth', 2);
xlabel('x');
ylabel('f(x)');
title('Sparsification Distribution');
grid on;






%% Lasso Regression via ChatGPT

% X: T-by-N matrix of time series (T samples, N nodes)
X = errts';
[T, N] = size(X);

% Store adjacency / coefficients
B = zeros(N, N);   % LASSO coefficients


spread = 0.045;
holder = ones(size(spread));
for k = 1:length(spread)

Lambda = spread(k);      % Regularization parameter (tune this)
    
for i = 1:N
    % The target node
    y = X(:, i);

    % All other nodes
    idx_other = setdiff(1:N, i);
    X_other = X(:, idx_other);

    % Fit LASSO (lasso returns coefficients for a range of lambdas)
    [beta_all, fitinfo] = lasso(X_other, y, 'Lambda', Lambda);

    % Extract coefficients for the chosen lambda
    beta = beta_all(:, 1);   % column 1 corresponds to Lambda

    % Store back into matrix B
    B(i, idx_other) = beta';
end

% Make symmetric adjacency (AND-rule)
connMatrix = (B ~= 0) & (B' ~= 0);

% Display the connectivity using a standard imagesc
figure()
imagesc(connMatrix)
xlabel('Parcel')
ylabel('Parcel')
title('Functional Connectivity Matrix')

% Calculate the BOLD predictions using just the resting state fMRI
predicted = ones(1,333) * connMatrix;
deltaPredicted = eField * connMatrix; 

% Calculate the partial linear correlation when controlling for variables
% in the resting state fMRI. Gives us a sense of how much of these results
% are due to the E-Field data
partialCorr = partialcorr(coeffs', deltaPredicted', predicted');
pairwiseR = corr(coeffs', deltaPredicted', Rows="Pairwise");
pairwiseRNoEF = corr(coeffs', predicted', Rows="Pairwise");

holder(k) = partialCorr; 

end

figure()
plot(coeffs, deltaPredicted, '.', MarkerSize=10)
xlabel('TMS BOLD Response')
ylabel('Predicted Response')
title(['Partial Coefficient ', num2str(partialCorr)])







%% Code for QNC Presentation 2
%
% Representation of how one might go about optimizing a LASSO regression
% approach, using ChatGPT


%% -------------------------------------------------------
%  Generate random data
%  X: predictor matrix (T x P)
%  y: outcome vector  (T x 1)
% --------------------------------------------------------
rng(1);                  % for reproducibility
T = 200;                 % number of samples
P = 20;                  % number of predictors

X = randn(T, P);         % Gaussian predictors
true_beta = [3; -2; 0; 0; 1; zeros(P-5,1)];
y = X * true_beta + 0.5*randn(T,1);   % linear model + noise


%% -------------------------------------------------------
%  Perform 10-fold cross-validated LASSO
% --------------------------------------------------------
K = 10;

[B, FitInfo] = lasso( ...
    X, y, ...
    'CV', K, ...                 % K-fold cross-validation
    'Standardize', true ...      % recommended
);

% B is P x L matrix of coefficients for all lambdas
% FitInfo contains CV results


%% -------------------------------------------------------
%  Identify the best lambda
% --------------------------------------------------------
idxLambdaMinMSE = FitInfo.IndexMinMSE;   % lambda minimizing cross-val MSE
bestLambda = FitInfo.Lambda(idxLambdaMinMSE);

fprintf('Optimal lambda (min MSE) = %.4f\n', bestLambda);


%% -------------------------------------------------------
%  Extract the beta coefficients for the optimal lambda
% --------------------------------------------------------
beta_optimal = B(:, idxLambdaMinMSE);
intercept = FitInfo.Intercept(idxLambdaMinMSE);

disp('Optimal LASSO coefficients:')
disp(beta_optimal)


%% -------------------------------------------------------
%  Optional: Plot cross-validation error curve
% --------------------------------------------------------
lassoPlot(B, FitInfo, 'PlotType','CV');
title('Cross-Validated MSE vs. Lambda');


%% -------------------------------------------------------
%  Optional: Fit final model using only optimal lambda
% --------------------------------------------------------
[B_final, FitInfo_final] = lasso(X, y, 'Lambda', bestLambda);

disp('Final model coefficients (same as beta_optimal):')
disp(B_final)
