clear;
close all;

% set of recorded test cases
emily_onsets = [3.111, 5.072, 7.097, 9.200, 9.989, 11.226, 12.896, 13.251, 15.234, 17.198, 19.050, 21.014, 22.292, 23.004, 25.016, 26.911, 28.601, 28.935];
emily_pitches = ["C2", "Ds2", "Gs2", "D2", "G1", "C2", "As1", "A1", "D2", "G1", "C2", "G1", "C2", "F1", "E1", "A1", "G1", "Fs1"];
corner_pocket_onsets = [1.887, 2.293, 2.724, 3.097, 3.506, 3.868, 4.261, 4.664, 5.053, 5.480, 5.892, 6.251, 6.653, 7.046, 7.427, 7.844, 8.216, 8.612, 9.017, 9.410, 9.794, 10.211, 10.635, 11.037, 11.464, 11.869, 12.247, 12.649, 13.063, 13.456, 13.822, 14.095, 14.212, 14.595, 14.991, 15.400, 15.808, 16.213, 16.576, 17.009, 17.426, 17.819, 18.200, 18.605, 18.992, 19.394, 19.793, 20.195, 20.613, 21.018, 21.402, 21.801, 22.206, 22.636, 23.020, 23.388, 23.790, 24.187, 24.573, 24.972, 25.359, 25.755, 26.154, 26.563, 26.977, 27.373, 27.776, 28.178, 28.561, 28.961, 29.341, 29.777, 30.164];
corner_pocket_pitches = ["Fs2", "F2", "Ds2", "Cs2", "C2", "Gs2", "As1", "C2", "Cs2", "F2", "Gs2", "A2", "As2", "As1", "F2", "D2", "Ds2", "F2", "Fs2", "G2", "Gs2", "A2", "As2", "C3", "Cs3", "As2", "Gs2", "Fs2", "F2", "C3", "As2", "E2", "E2", "Ds2", "F2", "Fs2", "Ds2", "Gs2", "Gs2", "C2", "Ds2", "Cs2", "Cs2", "B1", "B1", "As1", "F2", "D2", "As1", "Ds2", "Gs2", "Ds2", "Gs2", "Ds2", "Gs2", "Ds2", "Gs2", "Cs2", "Gs2", "F2", "Cs2", "Cs3", "B2", "Gs2", "G2", "Gs2", "Fs2", "F2", "Ds2", "Cs2", "D2", "Ds2", "F2"];
moons_harsh_onsets = [1.067, 2.592, 2.805, 4.045, 4.472, 6.219, 7.959, 9.675, 11.400, 13.140, 14.423, 14.871, 17.953, 18.179, 18.435, 20.219, 21.468, 21.920, 23.429, 23.651, 25.052, 25.509, 28.075, 28.989];
moons_harsh_pitches = ["C2", "C2"' "C2", "C2", "D2", "B1", "A1", "A2", "D2", "B1", "B1", "C2", "C2", "Cs2", "D2", "B1", "B1", "C2", "C2", "C2", "C2", "D2", "Cs2", "Fs1"];

real_onsets = {emily_onsets;corner_pocket_onsets;moons_harsh_onsets};
real_pitches = {emily_pitches;corner_pocket_pitches;moons_harsh_pitches};
real_names = ["emily.wav", "corner_pocket.wav", "moons_harsh.wav"];
    
note_names =    ["E1", "F1", "Fs1", "G1", "Gs1", "A1", "As1", "B1",... 
                "C2", "Cs2", "D2", "Ds2", "E2", "F2", "Fs2", "G2", "Gs2", "A2", "As2", "B2",...
                "C3", "Cs3", "D3", "Ds3", "E3", "F3", "Fs3", "G3", "Gs3", "A3", "As3", "B3",...
                "C4", "Cs4", "D4", "Ds4", "E4", "F4", "Fs4", "G4"];

note2midi_add = 48 - find(note_names == 'C2');

topdir = '/Users/rorybrooks/Desktop/4th_year/Jazzbox';
addpath(genpath(topdir));

fp_cumulative = 0;
fn_cumulative = 0;
tp_cumulative = 0;

f1_prec_recall_results = zeros(length(real_names), 3);

for n = 1:length(real_names)

    note_onsets_ground_truth = real_onsets{n}';
    these_pitches = real_pitches{n};
    note_pitches_ground_truth = zeros(length(note_onsets_ground_truth), 1);
    for i = 1:length(note_onsets_ground_truth)
        note_pitches_ground_truth(i) = find(note_names == these_pitches{i}, 1, "first");
    end
    note_pitches_ground_truth = note_pitches_ground_truth + note2midi_add;
    file_name = real_names(n);

    disp(file_name) 
    
    % may need to add directory to this
    this_file = file_name;
    [music_vec,fs_old]=readwav(this_file);
    
    % stereo handling
    if size(music_vec, 2) > 1
        music_vec = sum(music_vec, 2);
    end
    
    fs = 16000;
    NW = 256;
    overlap = 0.5;
    music_vec = v_resample(music_vec, fs, fs_old);
    
    %% compute onsets via function
    [onset_bool, offset_bool, onset_times, offset_times] = onset_detect_fn(music_vec, NW, overlap, fs);
    
    note_onsets_detected = onset_times;
    
    % initialise two storage arrays - one for the closest note onset at 
    note_onsets_closest = zeros(1, length(note_onsets_ground_truth));
    note_onsets_closest_diff = zeros(1, length(note_onsets_ground_truth));
    note_onsets_index = zeros(1, length(note_onsets_ground_truth));
    
    % so closest onset is the index of the closest note in ground truth for each entry in detected
    [onsets_diff, closest_onset] = min(abs(note_onsets_detected - note_onsets_ground_truth));
    for i = 1:length(note_onsets_ground_truth)
    
        this_note_options_ind = find(closest_onset == i);
        if isempty(this_note_options_ind)
            % therefore this note onset has not been caught
            note_onsets_closest(i) = nan;
            note_onsets_closest_diff(i) = nan;
            note_onsets_index(i) = nan;
        else
            [lat_diff, onset_index] = min(onsets_diff(this_note_options_ind));
    
            note_onsets_index(i) = this_note_options_ind(onset_index);
            note_onsets_closest(i) = note_onsets_detected(this_note_options_ind(onset_index));
            note_onsets_closest_diff(i) = onsets_diff(this_note_options_ind(onset_index));
        end
    end
    
    % note detection success = within 50ms
    notes_detect_success = (note_onsets_closest_diff <= 0.05);
    notes_detect_failure = (note_onsets_closest_diff > 0.05) | (isnan(note_onsets_closest_diff));
    notes_correctly_found_pc = 100 * (sum(notes_detect_success) / length(notes_detect_success));
    
    % correct detections - i.e. true positives
    correct_detections = max(length(find(notes_detect_success)), 0);
    tp_cumulative = tp_cumulative + correct_detections;
    tp = correct_detections;
    % over detections - i.e. false positives
    error_over_detections = max(length(note_onsets_detected) - length(note_onsets_ground_truth), 0);
    fp_cumulative = fp_cumulative + error_over_detections;
    fp = error_over_detections;
    % under detections - i.e. false negatives
    error_under_detections = max(length(find(notes_detect_failure)), 0);
    fn_cumulative = fn_cumulative + error_under_detections;
    fn = error_under_detections;
    
    %% now focus on the freq proc side of things:
    % turn frequency finding into a function, with input note onsets
    
    % run the pitch detection function - set the bayesian parameters here
    [autocorr_notes, note_estim_results] = wfl_pitch_detect(music_vec, fs, onset_bool, NW, overlap, 20, 5, 5, 0);
    
    % if needed, we can reduce the note_pitches_detected table to just those for correctly placed onsets
    valid_indices = find(~isnan(note_onsets_index) & notes_detect_success);
    valid_compar_max = note_estim_results(note_onsets_index(valid_indices), :);
    note_pitches_ground_truth_valid = note_pitches_ground_truth(valid_indices);
    
    % in midi pitch, 28 = E1
    % however, for our purposes, E1 = 1
    % so find midi pitch by adding 27
    note2midi_add = 48 - find(note_names == 'C2');
    
    note_pitches_detected_xcorr = autocorr_notes(note_onsets_index(valid_indices)) + note2midi_add;
    note_pitches_detected_wfl = valid_compar_max + note2midi_add;
    
    xcorr_accuracy = (note_pitches_detected_xcorr' == note_pitches_ground_truth_valid);
    
    
    % testing the waveform library is our main concern:
    
    wfl_first_detect_times = zeros(1, size(note_pitches_detected_wfl, 1));
    for i = 1:length(wfl_first_detect_times)
        this_note = note_pitches_ground_truth_valid(i);
    
        first_detected = find(note_pitches_detected_wfl(i, :) == this_note, 1, 'first');
        if isempty(first_detected)
            first_detected = nan;
        end
    
        wfl_first_detect_times(i) = first_detected;
        % this gives the first index where the wfl has correctly predicted the right note
    end
    
    
    wfl_error = note_pitches_detected_wfl - note_pitches_ground_truth_valid;
    wfl_error_all_notes = sum(wfl_error ~= 0, 1);
    
    wfl_error_note_times = zeros(1, size(wfl_error, 1));
    for i = 1:size(wfl_error, 1)
        wfl_error_note_times(i) = sum((wfl_error(i,:) ~= 0));
    end
 
    class_accuracy = sum(notes_detect_success == 1) / length(notes_detect_success);

    [f1_prec_recall_results(n, 1), f1_prec_recall_results(n, 2), f1_prec_recall_results(n, 3)] = proc_errors(tp, fp, fn);

end

[proc_f1, proc_precision, proc_recall] = proc_errors(tp_cumulative, fp_cumulative, fn_cumulative);

f1_prec_recall_results = round(f1_prec_recall_results, 4);

figure
% plot(wfl_error')
plot(wfl_error_note_times)
ylim([0, 1 + max(wfl_error_note_times)])