function wishart_filter(input,output)

addpath('./cifti','./wishart');

data=cifti_read(input);

DEMDT=1;              % Use 1 if demeaning and detrending (e.g. a timeseries) or -1 if not doing this (e.g. a PCA series)
VN=1;                 % Initial variance normalization dimensionality
Iterate=2;            % Iterate to convergence of dim estimate and variance normalization
NDist=2;              % Number of Wishart Filters to apply (for most single subject CIFTI grayordinates data 2 works well)
T=size(data.cdata,2); % Number of time points

% Wishart filter
filtered=icaDim(data.cdata,DEMDT,VN,Iterate,NDist);

% Demean and unit variance
filtered_data=detrend(filtered.data','constant');
filtered_data=filtered_data./repmat(std(filtered_data),T,1);
filtered_data=filtered_data';

% Save to file
data.cdata=filtered_data;
cifti_write(data, output);

end