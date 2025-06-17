function [onset_bool, offset_bool, onset_times, offset_times] = onset_detect_fn(music_vec, NW, overlap, fs)

    win = hann(NW);
    hop = NW * (1 - overlap);

    %% enframe and conduct 1 sided DFT for freq domain representation
    [frames,times,w] = v_enframe(music_vec,win,hop, 'f');
    frames = abs(frames);
    frames = frames';

    %% enframe using the same method, but this time keep it in the time domain
    sig_win = ones(1,NW);
    [sig_frames,sig_times,w] = v_enframe(music_vec,sig_win,hop);
    sig_frames = sig_frames';

    %% find mean of last avg_num frames. incr_freqs is based off of this
    avg_num = 5;
    frames_mean = zeros(size(frames));
    for i = 1:avg_num
        frames_mean = frames_mean + circshift(frames,i,2);
    end
    frames_mean = frames_mean ./ avg_num;
    
    incr_freqs = frames > (5 .* frames_mean);
    incr_freqs = sum(incr_freqs, 1);
    
    % implemet triangle numbers (multiply each by its position as a non-zero element per row)
    incr_freqs = (incr_freqs + (incr_freqs + 1)) ./ 2;
    
    % find mean of last avg_num incr_freqs value
    incr_avg = zeros(1, length(incr_freqs));
    for i = 1:avg_num
        incr_avg = incr_avg + circshift(incr_freqs,i,2);
    end
    incr_avg = incr_avg ./ avg_num; 
    
    incr_freqs = incr_freqs - incr_avg;
    incr_freqs(incr_freqs < 0) = 0;

    %% implement onsets    
    %% THRESHOLDING
    incr_freqs_threshold = 30;
    event_bool = incr_freqs > incr_freqs_threshold;
    
    %% signal energy vs prev average method
    sig_energy = sum(sig_frames .^ 2, 1);
    sig_avg_num = 10;
    sig_energy_diff = zeros(1, sig_avg_num);   
    for i = (sig_avg_num + 1):length(sig_energy)
        avg_last_num = sum(sig_energy((i - sig_avg_num):(i - 1))) / sig_avg_num;
        sig_energy_diff = [sig_energy_diff, (sig_energy(i) - avg_last_num)];
    end
    
    %% sig energy diff mad threshold
    sig_energy_diff_median = movmedian(sig_energy_diff,[5 0]);
    % circshift so the median is following but not including
    sig_energy_diff_median = circshift(sig_energy_diff_median, 1, 2);

    mad_k_energy = 2;
    mad_energy = movmedian(abs(sig_energy_diff - sig_energy_diff_median), [5 0]);
    mad_threshold_energy = sig_energy_diff_median + (mad_k_energy * mad_energy);
    
    %% go through an implement offset detection rules:
    % a) an offset must follow an onset (2 offsets cannot occur consecutively)
    % b) the energy difference of an offset must be negative
    %% it might make sense to rejig this using a measurement other than sig energy diff
    offset_bool = zeros(1, length(event_bool));
    onset_bool = zeros(1, length(event_bool));
    last_note_onset = false;

    %% integrate with sig energy diff
    for i = 1:length(offset_bool)
        % if there is an event at this index 
        if i > 1
            if event_bool(i) && ~onset_bool(i - 1)
                % if the event has negative sig energy diff and the last note was an onset
                if (sig_energy_diff(i) < mad_threshold_energy(i)) && last_note_onset
                    if(sig_energy_diff(i + 1) > mad_threshold_energy(i))
                        % it's just confused, this is an onset
                        onset_bool(i) = 1;
                        last_note_onset = true;
                    else
                        % then it's on offset
                        offset_bool(i) = 1;
                        last_note_onset = false;
                    end
                else
                    % otherwise it's an onset
                    onset_bool(i) = 1;
                    last_note_onset = true;
                end
            end
        else
            if event_bool(i)
                % if the event has negative sig energy diff and the last note was an onset
                if (sig_energy_diff(i) < mad_threshold_energy(i)) && last_note_onset
                    if(sig_energy_diff(i + 1) > mad_threshold_energy(i))
                        % it's just confused, this is an onset
                        onset_bool(i) = 1;
                        last_note_onset = true;
                    else
                        % then it's on offset
                        offset_bool(i) = 1;
                        last_note_onset = false;
                    end
                else
                    % otherwise it's an onset
                    onset_bool(i) = 1;
                    last_note_onset = true;
                end
            end
        end
    end
    
    onset_times = zeros(1, length(find(onset_bool)));
    offset_times = zeros(1, length(find(offset_bool)));
    j_times = 1;
    for i = 1:length(event_bool)
        if onset_bool(i)
            % times is given at the centre of the frame, but officially detection time is at the end of the frame.
            % also, times is given in samples
            % so... we need to add on half the frame length in samples
            onset_times(j_times) = (times(i) + (NW/2)) / fs;
            j_times = j_times + 1;
        end
    end
end