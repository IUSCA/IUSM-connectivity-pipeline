function vol_denoised = f_denoise_T1(vol,beta,patchradius,searchradius,rician,verbose)


%% default values
if nargin<6
    verbose = 0;
end
if nargin<5
    rician = 0;
end
if nargin<4
    searchradius = 1;
end
if nargin<3
    patchradius = 1;
end
if nargin<2
    beta = 1;
end

%% Normalize intensity range to [0,256]

map = isnan(vol(:));
vol(map) = 0;
map = isinf(vol(:));
vol(map) = 0;
mini = min(vol(:));
vol = vol - mini;
maxi = max(vol(:));
vol = vol*256/maxi;

%% Noise estimation
[hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(vol, rician, verbose);

%% Denoising
vol_denoised = MRIDenoisingONLM(vol, hfinal, beta, patchradius,  searchradius, rician, verbose);
map = vol_denoised<0;
vol_denoised(map)=0;

%% Restore Original intensity range
vol_denoised = (vol_denoised*maxi)/256;
vol_denoised = vol_denoised + mini;
