function [QCoutput]=bet_cycle_mricron(path2data)
%                       BET_CYCLE_MRICRON
%   A wrapper script that opens the T1_fov_denoised with
%   T1_brain_mask_filled overlaid for all subjects in the input directory
%   in sequence.
%
%   Uses mricron for visualization (path to mricron needs to be exported in
%   the user ~/.bashrc).
%
%   User feedback is requested once mricron is closed and is stored as 
%   output in the QCoutput cell along with the subject ID for review.
%
% Contributors:
%           Evgeny Chumin, Indiana University School of Medicine, 2018

%%
subjects = dir(path2data);
subjects(1:2)=[]; %remove '.' and '..'

sequence = {'T1'};

QCoutput=cell.empty(length(subjects),0);

for session=1:length(subjects)
    
    subject = subjects(session).name;
    QCoutput{session,1}=subject;
    
    T1Folder = [path2data '/' subject '/' sequence{1} '/'];
    fileIn = fullfile(T1Folder,'T1_fov_denoised.nii');
    maskIN = fullfile(T1Folder,'T1_brain_mask_filled.nii.gz');
    
    fprintf('\n Subject %s \n',subjects(session).name);
    %mricron T1_fov_denoised.nii -b 60 -o T1_brain_mask_filled.nii.gz -x
    sentence = sprintf('mricron %s -b 60 -o %s -x',fileIn, maskIN);
    [~,~] = system(sentence);
    reply = input('Does the brain mask have good coverage? \n','s');
    QCoutput{session,2}=reply;
    disp(QCoutput)
end 