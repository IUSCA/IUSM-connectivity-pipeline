function ts_outliers = f_detect_univariate_outliers(ts_data,scrubbing)

ts_outliers = false(size(ts_data));
for i=1:size(ts_data,1)
    ts = ts_data(i,:);
    thUp = prctile(ts(scrubbing),85)+(1.5*iqr(ts(scrubbing)));
    thLow = prctile(ts(scrubbing),15)-(1.5*iqr(ts(scrubbing)));
    ts_outliers(i,:) = (ts>thUp) | (ts<thLow);
end