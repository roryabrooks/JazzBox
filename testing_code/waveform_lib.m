classdef waveform_lib
    %WWAVEFORM_LIB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lib 
        cand_lib
        cand_lib_used = false;
        lib_status
        lib_size
        wf_length
        threshold
        freqs
        periods
        
    end
    
    methods
        function WFL = waveform_lib(note_freqs, ms100_numsamps, fs, thresh_factor)
            %WAVEFORM_LIB constructor

            if ~exist('thresh_factor','var')
                thresh_factor = 0.6;
            end

            % if isempty(num_entries)
            %     num_entries = 1;
            % end
            % for the minute, we'll just ignore the idea of a larger
            % waveform library

            WFL.lib_size = length(note_freqs);
            WFL.freqs = note_freqs;
            WFL.periods = double(1 ./ note_freqs);
            WFL.wf_length = ms100_numsamps;
            WFL.lib = zeros(ms100_numsamps, WFL.lib_size);
            WFL.lib_status = zeros(1,WFL.lib_size);

            dt = 1/fs;
            StopTime = (ms100_numsamps - 1) / fs; 
            t = (0:dt:StopTime)';

            for i = 1:WFL.lib_size
                this_sine = sin(2 * pi * t .* note_freqs(i));
                this_sine = this_sine ./ rms(this_sine);
                WFL.lib(:,i) = this_sine';
                WFL.lib_status(i) = -1;
            end

            WFL.threshold = floor(ms100_numsamps * thresh_factor);

        end

        function WFL = waveform_lib_update(WFL, update_waveform, update_note_index, wfl_detect_index, update_waveform_strength, update_rule) 
            %WAVEFORM_LIB_UPDATE update waveform library according to update rule

            length_diff = WFL.wf_length - length(update_waveform);
            if length_diff > 0
                update_waveform(end + 1:WFL.wf_length) = zeros(1,length_diff);
            elseif length_diff < 0
                update_waveform = update_waveform(1:WFL.wf_length);
            end

            if abs(rms(update_waveform) - 1) > 1e-6
                % disp("updating lib - changing rms")
                update_waveform = update_waveform ./ rms(update_waveform);
            end

            switch update_rule
        
                case 1 % rule 1 = naive update
                    WFL.lib(:,update_note_index) = update_waveform;
                    WFL.lib_status(update_note_index) = 1;
        
                case 2 % rule 2 = threshold update
                    if (update_waveform_strength > WFL.threshold) && (abs(update_note_index - wfl_detect_index) <= 1)
                        % disp("surpassed thresh! - updating wfl")
                        update_note_index = wfl_detect_index;
                        WFL.lib(:,update_note_index) = update_waveform;
                        WFL.lib_status(update_note_index) = 1;
                    else      
                        % disp("xcorr and wfl don't agree/ thresh not surpassed - no update")
                        % disp(update_note_index)
                        % disp(wfl_detect_index)
                    end
        
                case 3 % rule 3 = competitive update - replace worst
                    if WFL.cand_lib_used == false 
                        WFL.cand_lib = zeros(size(WFL.lib));
                        WFL.cand_lib_used = true;
                    end

                    % we first need to establish which of the library or
                    % candidate library entries correlate better: and thus
                    % which we will be replacing

                    if update_waveform_strength > WFL.threshold
                        % STILL SUBJECT TO THRESHOLD LIMIT!
               
                        %% perform cross correlation
                        % xcorr takes the second wave and shifts it from left to right over the first
                        compar_main = xcorr(update_waveform, WFL.lib(:,update_note_index));
                        compar_cand = xcorr(update_waveform, WFL.cand_lib(:,update_note_index));
    
                        max_compar_main = max(compar_main);
                        max_compar_cand = max(compar_cand);
    
                        if max_compar_cand > max_compar_main
                            % if the candidate is a better fit than the current lib entry
                            % place cand into lib, and place update as new cand
                            WFL.lib(:,update_note_index) = WFL.cand_lib(:,update_note_index);   
                        end
                        % if candidate is not a better fit, we just keep the current lib entry
                        % either way, place update into candidate
                        WFL.cand_lib(:,update_note_index) = update_waveform;
                    end

                % case 4 % rule 4 = competitive update - replace 2nd best
                % replace second best seems like a stupid idea for the
                % moment
        
                otherwise
                    % disp("invalid update rule given")
            end
        end

        function WFL = waveform_lib_interpolate(WFL)
            %WAVEFORM_LIB_INTERPOLATE interpolates the currently filled entries to the unfilled ones by stretching

            needs_update = find(WFL.lib_status == -1);
            filled_ins = find(WFL.lib_status ~= -1);

            if isempty(filled_ins) || isempty(needs_update)
                % can't be interpolated if everything is unfilled (only default sine waves)
                % OR if everything is filled
                % disp("error - library either completely full or empty - interp not possible")
                return 
            else
        
                num_empty = length(needs_update);
                num_filled = length(filled_ins);
        
                if num_filled < 1
                    % disp("cannot update - nothing filled in yet")
                    return
                elseif num_filled > 1
                    % interp1 needs at least 2 points. 
                    updaters = interp1(filled_ins, filled_ins, needs_update, 'nearest', 'extrap');
                else
                    % If we only have 1 then interp1 is useless
                    updaters = repelem(filled_ins, num_empty);
                end
        
                for i = 1:length(updaters)
        
                    this_period = WFL.periods(updaters(i));
                    desired_period = WFL.periods(needs_update(i));
                    [p, q] = rat(desired_period / this_period, 1e-6);
        
                    waveform_updater = resample(WFL.lib(:,updaters(i)), p, q);
        
                    length_diff = WFL.wf_length - length(waveform_updater);
        
                    if length_diff > 0  
                        % if updater is too short
                        % fill in remaining space with zeros
                        waveform_updater(end + 1:WFL.wf_length) = zeros(1,length_diff);
                    else 
                        % if updater is too long
                        waveform_updater = waveform_updater(1:WFL.wf_length);
                    end
        
                    WFL.lib(:,needs_update(i)) = waveform_updater;
                end
            end
        end
    end
end

