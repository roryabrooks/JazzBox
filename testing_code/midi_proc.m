clear;
close all;

function [f1, precision, recall] = proc_errors(tp, fp, fn)
    precision = tp / (tp + fp);
    recall = tp / (tp + fn);
    f1 = 2 * ((precision * recall) / (precision + recall));
end

addpath(genpath('/Users/rorybrooks/Desktop/4th_year/Jazzbox'));

note_matrix = readmidi('misty_bass_short.mid');
% this produces a note matrix, which has the form:
% ONSET(BEATS) | DURATION(BEATS) | MIDI channel | MIDI PITCH | VELOCITY | ONSET (SEC) | DURATION (SEC)  




%% note onset processing
note_onsets_ground_truth = note_matrix(:, 6);

note_onsets_detected = onset_times;
% for this_compar, each row is a note, and each column is an increasing window

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
% over detections - i.e. false positives
error_over_detections = max(length(note_onsets_detected) - length(note_onsets_ground_truth), 0);
% under detections - i.e. false negatives
error_under_detections = max(length(find(notes_detect_failure)), 0);





%% now focus on the freq proc side of things:
% turn frequency finding into a function, with input note onsets
note_pitches_ground_truth = note_matrix(:, 4);

note_names =    {'E1', 'F1', 'Fs1', 'G1', 'Gs1', 'A1', 'As1', 'B1',... 
                'C2', 'Cs2', 'D2', 'Ds2', 'E2', 'F2', 'Fs2', 'G2', 'Gs2', 'A2', 'As2', 'B2',...
                'C3', 'Cs3', 'D3', 'Ds3', 'E3', 'F3', 'Fs3', 'G3', 'Gs3', 'A3', 'As3', 'B3',...
                'C4', 'Cs4', 'D4', 'Ds4', 'E4', 'F4', 'Fs4', 'G4'};

% frequency processing:
[autocorr_notes, note_estim_results] = wfl_pitch_detect(music_vec, )

% if needed, we can reduce the note_pitches_detected table to just those for correctly placed onsets
valid_indices = find(~isnan(note_onsets_index) & notes_detect_success);
valid_compar_max = this_compar_max(note_onsets_index(valid_indices), :);
note_pitches_ground_truth_valid = note_pitches_ground_truth(valid_indices);

% in midi pitch, 28 = E1
% however, for our purposes, E1 = 1
% so find midi pitch by adding 27
note2midi_add = 48 - find(note_names == 'C2');

note_pitches_detected_xcorr = autocorr_notes(note_onsets_index(valid_indices));
note_pitches_detected_wfl = valid_compar_max + note2midi_add;

% testing the waveform library is our main concern:

wfl_first_detect_times = zeros(size(note_pitches_detected_wfl, 1));
for i = 1:length(wfl_first_detect_times)
    this_note = note_pitches_ground_truth_valid(i);

    disp(find(note_pitches_detected_wfl(i, :) == this_note, 1, 'first'))

    wfl_first_detect_times(i) = find(note_pitches_detected_wfl(i, :) == this_note, 1, 'first');
    % this gives the first index where the wfl has correctly predicted the right note
end


wfl_error = note_pitches_detected_wfl - note_pitches_ground_truth_valid;

figure
plot(wfl_error')


