function [results] = get_ndRMS_toShare(data,layout,IED,sample_freq,bad_channels,varargin)
%%% requirements
%add the latest version of FieldTrip toolbox
%(https://www.fieldtriptoolbox.org/) to your path.

%%% ------------ inputs:
% data: raw data [channels x time]
% layout: electrode location in the grid (e.g. [1,2,3;4,5,6;7,8,9;10,11,12]
% IED: inter-electrode distance
% sample_freq = sampling frequency
% bad_channels: bad channels selection
% optional:
% a matrix containing the start and stop events for each trial, where each
% row is a trial.

%%% ------------ outputs:
%freqs: the frequencies over which you computed the ndRMS
%ndRMS: the ndRMS
%IEDs: the inter-electrode distances the ndRMS was calculated for


plot_3d_plot = 1; %if you'd like to see a plot with ndRMS as a function of IED.

%% Create distance matrix
% careful, this approach only works for square or rectangular arrays. For
% other shapes (e.g. hexagonal), make sure to create the correct distance
% matrix in another way.

if size(varargin,2)>0
    if size(varargin{1},2) == 2
        trial_info = varargin{1};
        trial_flag = 1;
        fprintf('Found start stop events for %d trials...\n',size(varargin{1},1))
    else
        error('Trial information entered in the wrong dimensions.')
    end
else
    trial_flag = 0;
    fprintf('No trial information entered...\n')
end

[row,col] = deal(zeros(1,max(layout(:))));
for i = 1:max(layout(:)) %find index of each electrode
    [row(i), col(i)] = find(layout == i);
end

pair = [row; col]';
D = pdist(pair, 'euclidean');
Z = squareform(D); %create distance matrix
Z_adj = Z * IED;
Z_adj = round(Z_adj,3);
unique_distances = unique(Z_adj);
unique_distances(1) = []; %remove 0 distance

%% prepare analysis
chansel_corr = 1:size(data,1);
chansel_corr(bad_channels) = []; %to exlude channels from analysis.
g=sprintf('%d ', bad_channels);
fprintf('Bad channels %sare excluded from analysis\n', g)

% define fieldtrip hdr struct
hdr.Fs = sample_freq;
hdr.nchan = size(data,1);
hdr.nSamples = size(data,2);
hdr.nTrials = 1;
hdr.label = arrayfun(@num2str, 1:size(data,1),'UniformOutput', false)';

%cfg struct:
cfg = []; %empty again.
cfg.hdr                 = hdr;
cfg.trl                 = [1 length(data) 0 1];

% processing
trial = {};
time  = {};
trialinfo = 0;

begsample   = cfg.trl(1,1);
endsample   = cfg.trl(1,2);
trial{1}    = data(:,begsample:endsample);
begtime     = cfg.trl(1,3)/hdr.Fs;
endtime     = (cfg.trl(1,2)-cfg.trl(1,1)+cfg.trl(1,3))/hdr.Fs;
time{1}     = begtime:1/hdr.Fs:endtime;

p_data            = [];
p_data.hdr        = hdr;
p_data.fsample    = hdr.Fs;
p_data.sampleinfo = cfg.trl(:,1:2);
p_data.trialinfo  = trialinfo;
p_data.trial      = trial;
p_data.time       = time;
p_data.label      = hdr.label;
p_data.cfg        = cfg;

% Set channel selection and bad channels from subj struct
%first, bandstop filter.
cfg = [];
cfg.bsfilter   = 'yes';
cfg.bsfiltord  = 4;


ln_fr = 50;  %line noise frequency
cfg.bsfreq     = [((ln_fr:ln_fr:950)-1)' ((ln_fr:ln_fr:950)+1)'];

p_data         = ft_preprocessing(cfg,p_data);


%% unipolar processing
%1) re-referrencing for unipolars, if you want
disp('re-referencing data...')

cfg = [];
%cfg.refmethod     = 'avg';
cfg.reref       = 'yes';
cfg.refmethod   = 'median'; %common median re-referncing
cfg.refchannel  = chansel_corr; %rereference exculding bad channels.

p_data          = ft_preprocessing(cfg,p_data);
dat_input       = p_data.trial{1}';

% 2) extract frequency bins
fr_oi = [0 0; 1 4; 4 8;8 16; 16 32; 32 64; 64 128;128 256; 256 499;];
order = [nan; 3; 3; 4;4;4;4;4;4];

%set
data_filt_uni = nan(length(dat_input),size(data,1));

for fr = 1:size(fr_oi,1)

    if fr_oi(fr,1) == 0 && fr_oi(fr,2) == 0 %indicates unfiltered data. don't filter here.
        fprintf('Starting with unfiltered data...\n')
        data_filt_uni(:,chansel_corr) = dat_input(:,chansel_corr);
    else

        fprintf('Starting with %d-%d Hz...\n',fr_oi(fr,1),fr_oi(fr,2))

        disp('Filtering...')
        [b, a] = butter(order(fr), [fr_oi(fr,1)/(sample_freq/2) fr_oi(fr,2)/(sample_freq/2)], 'bandpass'); %try 2 or 3.
        for i = chansel_corr(1:end)
            %fprintf('filtering unipolar %d\n',i)
            %bandpass
            to_filt = dat_input(:,i);
            data_filt_uni(:,i) = filtfilt(b, a, to_filt);

        end
    end

    %%if necessary, extract trials
    if trial_flag == 1
        for i = 1:size(trial_info,1)
            task_state_nz(i,:,:) = data_filt_uni(trial_info(i,1):trial_info(i,2),:);
        end
    elseif trial_flag == 0
        task_state_nz(1,:,:) = data_filt_uni;
    end

    %zscore per trial.
    task_state = normalize(task_state_nz,2);

    disp('Unipolar processing finished...')

    %% create differentials

    disp('Starting bipolar processing...')

    %initialize differentials
    rms_diffs = nan(size(task_state,1),size(task_state,3),size(task_state,3));

    for tr = 1:size(task_state,1) %per trial
        diffs = nan(size(task_state,2),size(task_state,3),size(task_state,3));
        intermediate = squeeze(task_state(tr,:,:)); %bring down dimensions to [time x channels]
        for i = chansel_corr(1:end-1)
            ind = find(chansel_corr == i) + 1; %subtract only next channels (to avoid doubles and autocorrelations)

            for j = chansel_corr(ind:end)
                % fprintf('subtracting %d from %d\n',j,i)

                diffs(:,i,j) = intermediate(:,i) - intermediate(:,j);

            end
        end

        %rms over time. resulting dimensions: [trial x channel i x channel j]
        rms_diffs(tr,:,:) = squeeze(rms(diffs,1,'omitnan'));
    end


    %Median over trials)
    median_rms_diffs = squeeze(nanmedian(rms_diffs,1));

    %% group per distance
    disp('Grouping per distance...')
    min_number = 9;% minimal number of electrode pairs
    yy = 0;
    diff_p_distance = cell(size(unique_distances,1),1);
    for m = 1:size(unique_distances,1) %for each unique distance value
        diff_p_distance{m,:} = median_rms_diffs(Z_adj==unique_distances(m)); %find all values that correspond to a certain distance value
        num_obs = sum(~isnan(diff_p_distance{m,:}));
        if num_obs > min_number
            yy = yy + 1;
            diff_p_distance_mean(yy) = nanmean(diff_p_distance{m,:}); %#ok<*AGROW,*NANMEAN>
            IEDs(yy) = unique_distances(m);
        end
    end

    %% create output
    results(1,fr).freqs = [fr_oi(fr,1) fr_oi(fr,2)];
    results(1,fr).ndRMS = diff_p_distance_mean;
    results(1,fr).IEDs = IEDs;

end
%% plot 3D plot
if plot_3d_plot == 1
    figure('units','normalized','outerposition', [0 0 1 1]); hold on;

    Z = vertcat(results.ndRMS);
    X = vertcat(results(1).IEDs);
    X = repmat(X,[size(fr_oi,1) 1]);
    Y = meshgrid(1:size(fr_oi,1),X(1,:))';

    s = mesh(X,Y,Z);

    set(s,'LineWidth',2)
    shading interp
    xlabel('IED (mm)')
    zlabel('ndRMS (a.u.)')
    titles={'Unfiltered';'1-4';'4-8';'8-16';'16-32';'32-64';'64-128';'128-256';'256-499'};
    yticklabels(titles)

    view([-140,45])

    zlim([0.6 1.7])
    clim([0.6 1.7])

    set(gca,'FontSize',25,'Color','w')
    grid on;
end


end
