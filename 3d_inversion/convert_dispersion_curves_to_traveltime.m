clear;
clc;

corr_path = './acausal_correlations_dat/TT';
final_curves_dir = './disp_cur_TT/final_curves/all';
output_dir = './observed_times/TT';

% DO NOT EDIT BELOW THIS LINE
% --------------------------------
files = dir(final_curves_dir);
files = files(~ismember({files.name},{'.','..'}));

% loop over final inversion curves
for i = 1:length(files)
    % read the automatically picked curve since it contains the
    % coordinates of the involved stations
    fname_corr = erase(files(i).name,"CDisp.T.");
    fname_corr = erase(fname_corr, "GDisp.");
    fname_corr = fullfile(corr_path, fname_corr);

    X = importdata(fname_corr);
    
    % read station coordinates
    lon1 = X(1,1);
    lat1 = X(1,2);
    elv1 = X(1,3);

    lon2 = X(2,1);
    lat2 = X(2,2);
    elv2 = X(2,3);
    
    % compute station-station distance
    dis = distance(lat1, lon1, lat2, lon2);
    dis = deg2km(dis);
    
    % station elevation difference
    staelvdiff = abs(elv2 - elv1) / 1000;
    
    % correct station distance
    dis  = sqrt(dis*dis + staelvdiff*staelvdiff);
    
    % read final dispersion curve
    fname_final = fullfile(final_curves_dir, files(i).name);
    X = importdata(fname_final);
    
    % convert velocity to time
    X2 = X;
    X2(:,2) = dis ./ X(:,2);
    
    % write observed times
    fname_output = fullfile(output_dir, files(i).name);
    writematrix(X2, fname_output, 'Delimiter', '\t');
end

fprintf('done')
