% Demo for running linear descriminant analysis on example data and
% plotting ROC curves with confidence intervals and AUC

% B. Irving

% load an example dataset. The fisheriris dataset contains 3 labels
% ('virginica, versicolor, senosa) and features from thos species 
% http://www.mathworks.co.uk/products/demos/statistics/classdemo.html

load fisheriris
%loads species and meas

% Classification 1
% creating [0, 1] labels for classifying virginica or not virginica
labels = strcmp(species, 'virginica');

% Initialise the classifier class with the features and labels
S2=SimpleClassifier(meas, labels, 'virginica');
% Run leave one out cross validation
S2=S2.clas_LOOCV('linear'); 
% Display the confusion matrix of the LOOCV
S2.disp_conf();
% Generate a ROC curve
S2=S2.roc_curve_perf_pos();

% Classification 2
% creating [0, 1] labels for classifying versicolor or not versicolor
labels = strcmp(species, 'versicolor');

S3=SimpleClassifier(meas, labels, 'versicolor');
S3=S3.clas_LOOCV('linear'); 
S3.disp_conf();
S3=S3.roc_curve_perf_pos();

% Plot as many ROC objects as you want
plot(S2, S3)
