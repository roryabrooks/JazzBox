function [xcorr_notes, note_estim_results] = wfl_pitch_detect(music_vec, fs, onset_vector, NW, overlap, buffer_long, buffer_short, notes_considered, bayesian_mode)

    win = ones(1,NW); 
    hop = NW * (1 - overlap);
    [sig_frames,sig_times,w] = v_enframe(music_vec,win,hop);
    sig_frames = sig_frames';

    ms100_bound = ceil((((fs / 10) - NW) / (NW * (1 - overlap))) + 1); 
    ms100_numsamps = ((ms100_bound + 1) * (NW * (1 - overlap)));

    %% initialise note names
    note_names =    {'E1', 'F1', 'Fs1', 'G1', 'Gs1', 'A1', 'As1', 'B1',... 
                    'C2', 'Cs2', 'D2', 'Ds2', 'E2', 'F2', 'Fs2', 'G2', 'Gs2', 'A2', 'As2', 'B2',...
                    'C3', 'Cs3', 'D3', 'Ds3', 'E3', 'F3', 'Fs3', 'G3', 'Gs3', 'A3', 'As3', 'B3',...
                    'C4', 'Cs4', 'D4', 'Ds4', 'E4', 'F4', 'Fs4', 'G4'};
    
    note_names = convertCharsToStrings(note_names);

    %% initialise note frequencies and periods
    % note freqs in Hz - under 12TET and A = 440
    note_freqs_base =    [41.2034, 43.6535, 46.2493, 48.9994, 51.9131, 55.0000, 58.2705, 61.7354, 65.4064, 69.2957, 73.4162, 77.7817];   
    note_freqs = [];
    for i = 0:3
        note_freqs = [note_freqs, note_freqs_base .* (2 ^ i)];
    end
    note_freqs = note_freqs(1:length(note_names));
    note_periods = double(1 ./ note_freqs);

    %% initialise filling out waveform library with sine waves
    wf_lib = waveform_lib(note_freqs, ms100_numsamps, fs);

    %% initialise the note buffer
    note_hist = note_buffer(buffer_long, buffer_short, note_names);

    %% finding min possible lag for large scale autocorrelation
    % this = period of highest note possible / period of sampling freq
    % max possible lag, with the leeway of the distance between the longest and second longest periods
    min_poss_lag = floor(min(note_periods) * fs);
    max_poss_lag = ceil((max(note_periods) + (max(note_periods) - max(note_periods(note_periods~=max(note_periods))))) * fs);

    %% go through the signal and perform pitch ID
    % onset counter, used for indexing
    onset_counter = 1;
    this_compar = {};
    this_compar_max = [];
    
    % autocorr note detection storing
    autocorr_notes = [];
    autocorr_periods = [];
    autocorr_strengths = [];

    % initialise confidence store
    k_conf_inf = 0;
    bayesian_dist = zeros(1, length(note_names));

    % for every frame
    for i = 1:length(sig_frames)
    
        % if an onset is identified in this frame
        if onset_vector(i) 
    
            % reset comparison storage for this onset
            this_compar_times = cell(1, ms100_bound);
            this_time_max_note = zeros(1, ms100_bound);
            this_time_max_index = zeros(1, ms100_bound);
            
            % iterate through lengthening time window
            for t = 0:(ms100_bound - 1)
                if t == 0
                    analysis_win = sig_frames(:,i)';
                else
                    if i + t <= length(sig_frames)
                        analysis_add = sig_frames((hop+1):end, i + t)';
                        analysis_win = [analysis_win, analysis_add];
                    else
                        disp(["window overruns the signal! Stopped at window iteration: " t])
                        break
                    end
                end
    
                % signal length at this iteration
                % you could also just do sig_length = length(analysis_win);
                sig_length = NW + (t * NW * (1 - overlap));
                
                %% normalise analysis win by RMS! 
                this_analysis = analysis_win ./ rms(analysis_win);
    
                rms_this_analysis = rms(this_analysis);
                plot_rms = repmat(rms_this_analysis, length(analysis_win));
    
                % reset comparison storage for this note
                this_compar_notes = zeros(1,wf_lib.lib_size);
                this_compar_lags = zeros(1,wf_lib.lib_size);
    
                for w = 1:wf_lib.lib_size
                    %% extract relevant wave from waveform library and normalise by it's RMS
                    this_wave = wf_lib.lib(1:sig_length,w);
                    this_wave = (this_wave ./ rms(this_wave))';
                    
                    %% perform cross correlation
                    [compar, lags] = xcorr(this_analysis, this_wave);
    
                    [max_compar, max_ind] = max(compar);
                    max_lag = lags(max_ind);
    
                    this_compar_notes(w) = max_compar;
                    % for now we're ignoring the lag - it seems largely irrelevant in this context
                    this_compar_lags(w) = max_lag;
                end
    
                % normalise by the mean so that comparison is easier across scales (relative to 1, which will be mean)
                this_compar_notes = this_compar_notes ./ mean(this_compar_notes);

                conf_wfl = 1 - (1 / (max(this_compar_notes)/mean(this_compar_notes)));
                % bayesian_dist
                % k_conf_inf

                switch bayesian_mode
                    case 0 
                        % 0 = don't use
                    case 1 
                        % 1 = use naive - 0.5
                        this_compar_notes = (0.5 .* this_compar_notes) + (0.5 .* bayesian_dist);
                    case 2
                        % 2 = use kalman gain
                        this_compar_notes = (conf_wfl .* this_compar_notes) + (k_conf_inf .* bayesian_dist);
                    case 3
                        % 3 = use 0.1
                        this_compar_notes = (0.9 .* this_compar_notes) + (0.1 .* bayesian_dist);
                    otherwise
                        disp("bayesian mode invalid!")
                end
    
                this_compar_times{t + 1} = this_compar_notes;
                [this_time_max_note(t + 1), this_time_max_index(t + 1)] = max(this_compar_notes);
                                 
            end 
    
            %% storage of results
            this_compar = [this_compar; this_compar_times];
            this_compar_max = [this_compar_max; this_time_max_index];

            %% autocorrelation step
            % when the window is at it's largest, perform autocorrelation to detect true note
            [true_note, true_lags] = xcorr(this_analysis);

            % figure
            % plot(true_lags, true_note)
            % pause
            % close all;
               
            true_note_deriv_1 = true_note - circshift(true_note, 1);
            % choose turning points to be only maximums
            tp_locs = ((true_note_deriv_1 < 0) & (circshift(true_note_deriv_1, 1) > 0));
            tp_magns = true_note(tp_locs);
            tp_lags = true_lags(tp_locs);
    
            % calculate and remove the zero peak (anything below min possible lag) and find the true lag
            tp_invalid_ind = tp_lags < min_poss_lag;
            tp_invalid_ind = tp_invalid_ind | (tp_lags > max_poss_lag);
            tp_magns(tp_invalid_ind) = [];
            tp_lags(tp_invalid_ind) = [];
    
            [this_true_magn, this_true_lag_ind] = max(tp_magns);
            this_true_lag = abs(tp_lags(this_true_lag_ind)); % abs, since lag could be negative     
            autocorr_strengths = [autocorr_strengths, this_true_magn];
    
            % find the period/note - closest in set to our autocorr calculated value
            this_period = this_true_lag / fs;
            [~, this_period_ind] = min(abs(note_periods - this_period));
            this_period_closest = note_periods(this_period_ind);
            this_note = this_period_ind;
            this_note_name = note_names(this_period_ind);
    
            % update using xcorr strength threshold. This whole process could be encompassed in the method
            % disp(["onset counter = " onset_counter])
            wf_lib = wf_lib.waveform_lib_update(this_analysis,this_period_ind, this_time_max_index(end), this_true_magn,2);
            wf_lib = wf_lib.waveform_lib_interpolate();
    
            autocorr_notes = [autocorr_notes, this_note];
            autocorr_periods = [autocorr_periods, this_period_ind];

            % update note hist and assess confidence: this ranges from 0 to 1
            note_hist = note_hist.add(char(this_note_name));
            [note_hist, k_conf_inf] = note_hist.assess_confidence();

            % disp(onset_counter)

            if onset_counter > notes_considered
                bayesian_dist = note_hist.bayesian_inference(notes_considered, 2, 2);
            end
            % end of autocorrelation step
    
            onset_counter = onset_counter + 1;
        end
    end

    xcorr_notes = autocorr_notes;
    note_estim_results = this_compar_max;
end