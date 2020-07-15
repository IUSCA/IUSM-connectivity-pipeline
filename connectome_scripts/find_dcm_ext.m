function [dcm_ext]=find_dcm_ext(path2dcmdir)

% accepted dicom formats
f_ext = {'dcm','IMA'};

f_ext_out=nan(length(f_ext),1);
for ff=1:length(f_ext)
    list_dicoms = dir(fullfile(path2dcmdir,sprintf('*.%s',f_ext{ff})));
    if size(list_dicoms,1) > 0
        f_ext_out(ff,1) = 1;
    else
        f_ext_out(ff,1) = 0;
    end
end

% ensure only one extension is found
if sum(f_ext_out) == 1
    dcm_ext = f_ext{logical(f_ext_out)};
% if none or multiple found set to empty
else
    dcm_ext = [];
end