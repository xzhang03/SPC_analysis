function Areas = spcAUC(data)
% spcAUC performs area-under-curve analysis for data.
% Areas = spcAUC(data)

% Get dim
npts = size(data, 2);

% Slope from first to last point
slp = (data(:,end) - data(:,1)) / (npts - 1);

% Expected intermediate points
data_model = slp * (0 : (npts-1)) + data(:,1) * ones(1,npts);

% Deviation from model
Dev = data - data_model;

% Areas (sum of deviations)
Areas = sum(Dev,2);

end