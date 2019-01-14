function f_run_LPCA(DWI_fileIn,DWI_fileOut,rician,numThreads)
    
    beta=1;
    verbose = 0;
    
    if nargin<4
        numThreads = 6;
    end
    
    data = MRIread(DWI_fileIn);
    ima = data.vol;
    %s=size(ima);
    ima = double(ima);

    % fixed range
    map = isnan(ima(:));
    ima(map) = 0;
    map = isinf(ima(:));
    ima(map) = 0;
    mini = min(ima(:));
    ima = ima - mini;
    maxi=max(ima(:));
    ima=ima*255/maxi;          
    disp('before')
    size(ima)
    DWIdenoised = DWIDenoisingLPCA(ima,beta,rician,numThreads,verbose);
    disp('after');                  
    % Original intensity range
    DWIdenoised=DWIdenoised*maxi/255 + mini;
    map = find(DWIdenoised<0);
    DWIdenoised(map) = 0;
    map = isnan(DWIdenoised);
    DWIdenoised(map) = 0;
    data.vol = DWIdenoised;        
    MRIwrite(data,DWI_fileOut);
    
    

