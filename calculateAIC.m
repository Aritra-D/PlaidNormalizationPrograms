function aic = calculateAIC(SS,N,K)

%AIC Akaike's Information Criterion for curve-fitting model comparison
%   AIC(SS,N,K) where SS is the sum of squared residuals (as returned by e.g. lsqcurvefit)
%   N is the number of data points and K is the number of coefficients. 
%   Computes and returns the corrected AIC score for the model.
%   
%   The model with the lowest AIC value is the best fit.
%
%   References: (1) Motulsky, H. and Christopoulos, A. (2003) Fitting models to biological data using 
%   linear and nonlinear regression. A practical guide to curve fitting. Graphpad Software Inc., San Diego CA.
%   www.graphpad.com
%
%   (2) Akaike, H. (1974) "A new look at the statistical model identification." IEEE Transactions on Automatic Control, AC-19, 716-723
%
%   (3) Hurvich, C. M., and Tsai, C-L. (1989). Regression and time series model selection in small samples. Biometrika, 76, 297-307.
%   [the AIC correction]   
%
%   NOTE: this computes AIC from sum-of-squares (SS), and thus uses SS as
%   an esimator for the maximum likelihood; when actually fitting
%   distributions using MLE, then use AICL instead!
%
%   Mark Humphries 11/10/2004, ModelDB Accession: 128818

%  Humphries MD, Lepora N, Wood R, Gurney K (2009) 
% Capturing dopaminergic modulation and bimodal membrane behaviour of 
% striatal medium spiny neurons in accurate, reduced models. 
% Front Comput Neurosci 3:26

K  = K + 1; % additional degree-of-freedom is SS.

% raw AIC
aic = N .* log(SS./N) + 2 .* K;

% % apply correction in case N close to K
aic = aic + (2.*K.*(K+1)) ./ (N - K - 1);