function features  = get_features_json(json_file, sliceTiming_output, normalize)

% Obatin a struct of image features from a stored json files
% Will write out the slice timing table (n by 1) as a txt file. 

    if nargin < 3
        normalize = true;
    end
    
    if nargin < 2
        % save the slicetiming in the same folder
        sliceTiming_output = true;
    end
    
    features = jsondecode(fileread(json_file));
    
    % Normalize the sliceT if needed
    if normalize == true
        sliceT = features.SliceTiming ./ features.RepetitionTime - 0.5;
    else
        sliceT = features.SliceTiming;
    end 
    
    % Write out the sliceTiming into a text file
    parent_dir = strsplit(json_file, '/');
    parent_dir = strjoin(parent_dir(1:length(parent_dir)-1),'/');
    
    if sliceTiming_output
        writetable(table(sliceT), fullfile(parent_dir,'slicetimes_frac.txt'), ...
        'WriteVariableNames', false);
    end
 

end
