clear;
close all;

topdir = '/Users/rorybrooks/Desktop/4th_year/Jazzbox';
addpath(genpath(topdir));

% define which test set we want to use
synthdir = fullfile(topdir, 'sound_files/synthetic/jazz_set', '**');

files = dir(fullfile(synthdir, '*.wav'));
file_names = {files.name};
midis = dir(fullfile(synthdir, '*.mid'));
midi_names = {midis.name};

fp_cumulative = 0;
fn_cumulative = 0;
tp_cumulative = 0;
wfl_error_big = zeros(length(midi_names), 12);

f1_prec_recall_results = zeros(length(midi_names), 3);

for n = 1:length(midi_names)

    midi_name = string(midi_names(n));
    file_name = string(file_names(n));

    disp(file_name)
    
    note_matrix = readmidi(midi_name);
    % this produces a note matrix, which has the form:
    % ONSET(BEATS) | DURATION(BEATS) | MIDI channel | MIDI PITCH | VELOCITY | ONSET (SEC) | DURATION (SEC)  
    
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
    
    %% note onset processing
    note_onsets_ground_truth = note_matrix(:, 6);
    
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
    note_pitches_ground_truth = note_matrix(:, 4);
    
    note_names =    ["E1", "F1", "Fs1", "G1", "Gs1", "A1", "As1", "B1",... 
                    "C2", "Cs2", "D2", "Ds2", "E2", "F2", "Fs2", "G2", "Gs2", "A2", "As2", "B2",...
                    "C3", "Cs3", "D3", "Ds3", "E3", "F3", "Fs3", "G3", "Gs3", "A3", "As3", "B3",...
                    "C4", "Cs4", "D4", "Ds4", "E4", "F4", "Fs4", "G4"];
    
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

    wfl_avg = mean(abs(wfl_error), 1);

    wfl_error_big(n, :) = wfl_avg;
 
    class_accuracy = sum(notes_detect_success == 1) / length(notes_detect_success);

    [f1_prec_recall_results(n, 1), f1_prec_recall_results(n, 2), f1_prec_recall_results(n, 3)] = proc_errors(tp, fp, fn);



end

[proc_f1, proc_precision, proc_recall] = proc_errors(tp_cumulative, fp_cumulative, fn_cumulative);

f1_prec_recall_results = round(f1_prec_recall_results, 4);

% figure
% plot(wfl_error')
% plot(wfl_error_note_times)
% ylim([0, 1 + max(wfl_error_note_times)])

figure
plot(wfl_error_big')
title("Note Errors Present over Subsequent Frames of a Note in Jazz Test Set")
xlabel("Frame Index")
ylabel("Note Error (Semitones)")
% plot(1:(1/fs):music_vec)
% hold on
% plot(onset_times, ones(1, length(onset_times)))
% plot(wfl_error_note_times)
ylim([(min(wfl_error(:)) - 1), (1 + max(wfl_error(:)))])
xlim([1 size(wfl_error, 2)])

print('-depsc', 'filename.eps');