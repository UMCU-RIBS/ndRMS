%example script to test ndRMS function

%generate random data with mean 0 and std 1. 
data = randn(64, 200000);

%set required variables
sample_freq = 2000;
bad_channels = [1,2,3,4];
layout = [1:16;17:32;33:48;49:64]';
IED = 3; %mm

%optional: set trial info 
start_tr = 1:3000:175000;
stop_tr = start_tr + 2000;
trial_info = [start_tr;stop_tr]';

%run function
results = get_ndRMS(data,layout,IED,sample_freq,bad_channels,trial_info); 