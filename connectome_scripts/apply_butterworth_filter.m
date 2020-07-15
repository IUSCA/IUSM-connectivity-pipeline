function [tsf] = apply_butterworth_filter(ts,order,fmin,fmax)

[b,a] = butter(order,[fmin,fmax]);
tsf = filtfilt(b,a,ts);