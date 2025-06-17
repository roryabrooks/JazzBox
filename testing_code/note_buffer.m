classdef note_buffer
    %NOTE_BUFFER buffer for note history and key detection

    properties
        buffer_long
        buffer_short
        % these are the actual note histograms, counting the number of occurences of each note, positioned at C = 1 -> B = 12
        note_count_long
        note_count_short
        max_size_long
        max_size_short
        index_long
        index_short

        valid_notes = ["C", "Cs", "D", "Ds", "E", "F", "Fs", "G", "Gs", "A", "As", "B"]
        note_names
        maj_scale = [1, 3, 5, 6, 8, 10, 12]
        min_scale = [1, 3, 4, 6, 8, 9, 11]
        K_ind = [4, 11, 6, 1, 8, 3, 10, 5, 12, 7, 2, 9]
        KK_maj_prof = [6.35, 2.23, 3.48, 2.33, 4.38, 4.09, 2.52, 5.19, 2.39, 3.66, 2.29, 2.88]
        KK_min_prof = [6.33, 2.68, 3.52, 5.38, 2.60, 3.53, 2.54, 4.75, 3.98, 2.69, 3.34, 3.17]
        
        this_tonic_long
        this_is_maj_long
        this_tonic_short
        this_is_maj_short
        this_scale

        % note_movements records the movement type between the last N notes
        %  movement_dirs records whether the movement is up, down, or stationary for each
        note_movements
        index_movements
        movement_dirs
        max_size_movements
        last_note

        % AS_maj_prof = [7.0, 0.2, 2.3, 0.5, 2.7, 5.0, 0.4, 5.1, 0.4, 2.5, 0.3, 1.0];
        % AS_min_prof = [7.0, 0.6, 2.4, 5.1, 0.3, 4.8, 0.5, 4.3, 1.0, 0.7, 2.6, 0.5];
    end

    methods
        function obj = note_buffer(N_long, N_short, note_names)
            obj.max_size_long = N_long;
            obj.buffer_long = zeros(1, N_long);
            obj.note_count_long = zeros(1,12);
            obj.index_long = 0;

            obj.max_size_short = N_short;
            obj.buffer_short = zeros(1, N_short);
            obj.note_count_short = zeros(1, 12);
            obj.index_short = 0;

            obj.note_movements = NaN(1, N_long);
            obj.movement_dirs = zeros(1, N_long);
            obj.max_size_movements = N_long;
            obj.index_movements = 0;

            obj.note_names = note_names;
            obj.last_note = 0;

            obj.this_scale = zeros(1, 7);
        end

        function obj = add(obj, val)
            if isstring(val) | ischar(val)
                % if val is given as a note name
                find_str = val(1:end - 1);
                val_add = find(obj.valid_notes == find_str);
                this_note = find(obj.note_names == val);
                if isempty(this_note)
                    disp("Error - Note not in note names:")
                    disp(val)
                    return
                end
            else
                % otherwise, we can assume given as note index
                val_add = mod(val + 3, 12) + 1;
                this_note = val; 
            end

            % increment index
            obj.index_long = mod(obj.index_long, obj.max_size_long) + 1;
            % find val at this incremented index, and decrease its count in note histogram
            val_remove_long = obj.buffer_long(obj.index_long);
            if val_remove_long ~= 0
                obj.note_count_long(val_remove_long) = obj.note_count_long(val_remove_long) - 1;
            end
            % update note buffer and increment note histogram
            obj.buffer_long(obj.index_long) = val_add;
            obj.note_count_long(val_add) = obj.note_count_long(val_add) + 1;

            % do the same for the short list
            obj.index_short = mod(obj.index_short, obj.max_size_short) + 1;
            val_remove_short = obj.buffer_short(obj.index_short);
            if val_remove_short ~= 0
                obj.note_count_short(val_remove_short) = obj.note_count_short(val_remove_short) - 1;
            end
            obj.buffer_short(obj.index_short) = val_add;
            obj.note_count_short(val_add) = obj.note_count_short(val_add) + 1;

            % if the program has identified the key and scale
            if sum(obj.this_scale) ~= 0
                % find the note movement type from the last note to this one
                if obj.last_note ~= 0
                    obj.index_movements = mod(obj.index_movements, obj.max_size_movements) + 1;

                    note_movement = this_note - obj.last_note;
                    if note_movement > 0
                        note_dir = 1;
                    else
                        note_dir = max(note_movement, -1); 
                    end
                    
                    last_note_scale = find(obj.this_scale == (mod(obj.last_note + 3, 12) + 1));
                    this_note_scale = find(obj.this_scale == (mod(this_note + 3, 12) + 1));

                    if ~isempty(last_note_scale) & ~isempty(this_note_scale)
                    % in the case that both notes exist in the scale
                        if this_note_scale >= last_note_scale
                            scale_movement = this_note_scale - last_note_scale;
                        else
                            scale_movement = (this_note_scale + 7) - last_note_scale;
                        end
                        switch scale_movement
                        % scale is based on the indices of the scale list: [1, 2, 3, 4, 5, 6, 7]
                            case 0
                            % if they are the same scale degree
                                if obj.last_note ~= this_note
                                % but not the same note, we are moving in octaves
                                    movement_type = 8;
                                else
                                % else they are the same note
                                    movement_type = 0;
                                end
                            case 4
                            % movement type is in fifths
                                movement_type = 5;
                            case 3
                            % movement type is in fourths (combine with the fifths logic)
                                movement_type = 4;
                            case 2
                            % movement type is in thirds
                                movement_type = 3;
                            case 1
                            % movement type is by scale (2 semitones, so denoted 2)
                                movement_type = 2;
                            otherwise
                            % there is no recognised relationship between these scale tones
                                movement_type = -1;
                        end
                    elseif abs(note_movement) == 1 
                    % in the case that there is chromatic movement
                        movement_type = 1;
                    else
                    % so the two notes aren't both in the scale, and aren't chromatic
                    % therefore, movement type is undefined
                        movement_type = -1;   
                    end
            
                    % now add to the note_movements list
                    obj.note_movements(obj.index_movements) = movement_type;
                    obj.movement_dirs(obj.index_movements) = note_dir;

                end
            else
                disp("Error - no scale yet defined!")
            end
            obj.last_note = this_note;
        end

        function out = return_buffer_long(obj)
            % Return buffer in order: oldest to newest
            idx = mod((obj.index_long:obj.index_long + obj.max_size_long - 1), obj.max_size_long) + 1;
            out = obj.buffer_long(idx);
        end

        function out = return_buffer_short(obj)
            % Return buffer in order: oldest to newest
            idx = mod((obj.index_short:obj.index_short + obj.max_size_short - 1), obj.max_size_short) + 1;
            out = obj.buffer_short(idx);
        end

        function obj = return_key_SF_long(obj)
            % Return the key as detected by the signature of fifths method
            K = obj.note_count_long(obj.K_ind);
            K = K ./ max(K);

            axis_finder = ones(1, 6);
            K(end + 1: end + 5) = K(1:5);

            main_axis = zeros(1,12);

            for i = 1:12
                main_axis(i) = sum(axis_finder .* K(i:i+5));
            end
            % this lists the key likelihood in this order:
            % Ds, As, F, C, G, D, A, E, B, Fs, Cs, Gs
            % [4, 11, 6, 1, 8, 3, 10, 5, 12, 7, 2, 9]
            % ["C", "Cs", "D", "Ds", "E", "F", "Fs", "G", "Gs", "A", "As", "B"]

            % note - max notes is in fifths order, need to put in scale order
            max_notes = find(main_axis == max(main_axis));
            max_notes = obj.K_ind(max_notes);
            count_ext = [obj.note_count_long, obj.note_count_long];
            max_maj_corr = 0;
            max_min_corr = 0;
            max_maj_note = 0;
            max_min_note = 0;

            % disp(max_notes)
            
            for i = max_notes
                this_note_count = count_ext(i:i+11);
                maj_prof_corr = sum(obj.KK_maj_prof .* this_note_count);
                if maj_prof_corr > max_maj_corr
                    max_maj_corr = maj_prof_corr;
                    max_maj_note = i;
                end
                min_prof_corr = sum(obj.KK_min_prof .* this_note_count);
                if min_prof_corr > max_min_corr
                    max_min_corr = min_prof_corr;
                    max_min_note = i;
                end
            end

            if max_maj_corr >= max_min_corr
                obj.this_tonic_long = max_maj_note;
                obj.this_is_maj_long = true;
            else
                obj.this_tonic_long = max_min_note;
                obj.this_is_maj_long = false;
            end
        end   

        function obj = return_key_SF_short(obj)
            % Return the key as detected by the signature of fifths method
            K = obj.note_count_short(obj.K_ind);
            K = K ./ max(K);

            axis_finder = ones(1, 6);
            K(end + 1: end + 5) = K(1:5);

            main_axis = zeros(1,12);

            for i = 1:12
                main_axis(i) = sum(axis_finder .* K(i:i+5));
            end
            % this lists the key likelihood in this order:
            % Ds, As, F, C, G, D, A, E, B, Fs, Cs, Gs
            % [4, 11, 6, 1, 8, 3, 10, 5, 12, 7, 2, 9]
            % ["C", "Cs", "D", "Ds", "E", "F", "Fs", "G", "Gs", "A", "As", "B"]

            % note - max notes is in fifths order, need to put in scale order
            max_notes = find(main_axis == max(main_axis));
            max_notes = obj.K_ind(max_notes);
            count_ext = [obj.note_count_short, obj.note_count_short];
            max_maj_corr = 0;
            max_min_corr = 0;
            max_maj_note = 0;
            max_min_note = 0;

            % disp(max_notes)
            
            for i = max_notes
                this_note_count = count_ext(i:i+11);
                maj_prof_corr = sum(obj.KK_maj_prof .* this_note_count);
                if maj_prof_corr > max_maj_corr
                    max_maj_corr = maj_prof_corr;
                    max_maj_note = i;
                end
                min_prof_corr = sum(obj.KK_min_prof .* this_note_count);
                if min_prof_corr > max_min_corr
                    max_min_corr = min_prof_corr;
                    max_min_note = i;
                end
            end

            if max_maj_corr >= max_min_corr
                obj.this_tonic_short = max_maj_note;
                obj.this_is_maj_short = true;
            else
                obj.this_tonic_short = max_min_note;
                obj.this_is_maj_short = false;
            end
        end 

        function obj = get_scale(obj, type)
        % return the scale of the detected key, based on either the long history of short history (type = "long" or "short")

            if type == "long"
                this_tonic = obj.this_tonic_long;
                is_maj = obj.this_is_maj_long;
            elseif type == "short"
                this_tonic = obj.this_tonic_short;
                is_maj = obj.this_is_maj_short;
            else
                % disp("type must be given as (long) or (short)")
            end


            if is_maj
                scale = mod(obj.maj_scale + (this_tonic - 1), 12);
            else
                scale = mod(obj.min_scale + (this_tonic - 1), 12);            
            end
            scale(scale == 0) = 12;
            obj.this_scale = scale;
        end

        function [obj, k_conf] = assess_confidence(obj)

            %% assess our confidence in the key at this moment
            obj = obj.return_key_SF_short();
            obj = obj.return_key_SF_long();

            % update scale, just good practise
            obj = obj.get_scale("long");

            K_ind_dbl = [obj.K_ind, obj.K_ind];

            tonic_loc_short = find(K_ind_dbl == obj.this_tonic_short);
            tonic_loc_long = find(K_ind_dbl == obj.this_tonic_long);

            % Compute all pairwise absolute differences
            dist_matrix = abs(tonic_loc_long' - tonic_loc_short);  % A' is column, B is row → creates a matrix
            
            % Find the minimum value
            key_dist = min(dist_matrix(:));
            % 6 corresponds to a tritone, i.e. as far away as possible
            key_conf = 1 - (key_dist / 6);

            %% assess how much the note movements fit into our categorisations

            not_classifiable = sum(obj.note_movements == -1);
            movement_uncert = (not_classifiable / sum(~isnan(obj.note_movements)));
            movement_conf = 1 - movement_uncert;

            % independent probs multiply, we're assuming key conf and movement conf are indep
            k_conf = key_conf * movement_conf;
            
        end

        function out_dist = bayesian_inference(obj, n_notes, delay_falloff, oct_falloff)
        % return the output distribution of notes as defined by our probabalisitc model

            if(n_notes < 1)
                disp("Error - number of notes used for inference must be at least 1")
                return
            elseif n_notes > obj.max_size_long
                disp("Error - too many notes used for inference, greater than note buffer length")
                return
            elseif sum(~isnan(obj.note_movements)) < n_notes
                disp("Error - Not enough notes in movement history yet!")
                return
            end
    
            n_notes_indices = nan(1, n_notes);
            notes_split = obj.index_movements - n_notes;
            
            if notes_split >= 0
                n_notes_indices = ((obj.index_movements + 1) - n_notes:obj.index_movements);
            else
                n_notes_indices((end + 1) - obj.index_movements:end) = 1:obj.index_movements;
                n_notes_indices(1:end - obj.index_movements) = (((obj.max_size_long + 1) - abs(notes_split)):obj.max_size_long);
            end
    
            movement_hist = obj.note_movements(n_notes_indices);
            dir_hist = obj.movement_dirs(n_notes_indices);

            % disp("bayesian")
            % disp(movement_hist)
            % disp(n_notes_indices)

            if obj.this_is_maj_long
                diatonic_dist = obj.KK_maj_prof;
            else
                diatonic_dist = obj.KK_min_prof;
            end

            this_note = obj.buffer_long(obj.index_long);
            greater_than_tonic = (this_note + 1) - obj.this_tonic_long;
            if greater_than_tonic > 0
                this_degree = greater_than_tonic;
            else
                this_degree = greater_than_tonic + 12;
            end

            % replicate the key profile on either side of this degree:
            doubleside_dist = [diatonic_dist(this_degree:end), diatonic_dist(1:(this_degree - 1)), diatonic_dist(this_degree:end), diatonic_dist(1:this_degree)];
            % doubleside_dist runs from 1 to 25, and is centred with this note at 13
            % profile scaling sets the average magnitude of the doubleside dist
            doubleside_dist = (doubleside_dist ./ mean(doubleside_dist));
            last_plot_dist = nan(1,length(doubleside_dist));
            

            % calculate delay scaling according to delay falloff param
            delay_scaling = 1 ./ (delay_falloff .^ ((n_notes - 1):-1:0));

            up_degree = sum((dir_hist >= 0) .* delay_scaling);
            down_degree = sum((dir_hist <= 0) .* delay_scaling);
            % stat_degree = sum((dir_hist == 0) .* delay_scaling);

            degree_sum = up_degree + down_degree;
            up_degree = up_degree / degree_sum;
            down_degree = down_degree / degree_sum;
            % stat_degree = stat_degree / degree_sum;

            for i = n_notes:-1:1
                this_movement = movement_hist(i);

                % adjust the scale to ensure the tone numbers are always increasing
                this_scale_ind = find(obj.this_scale == (mod(this_note + 3, 12) + 1));
                this_scale_adjust_ind = obj.this_scale < [obj.this_scale(end), obj.this_scale(1:end - 1)];
                scale_adjust = obj.this_scale + 12 * this_scale_adjust_ind;

                % disp(obj.this_scale)
                % disp(scale_adjust)

                switch this_movement
                    case 8
                    % octave move
                        up_shift = 12;
                        down_shift = 12;
                    case 5
                    % movement by a fifth -> up by a fourth or down by a fifth next
                        up_shift = 5;
                        down_shift = 7;
                    case 4
                    % movement by a fourth -> up by a fourth or down by a fifth next
                        up_shift = 5;
                        down_shift = 7;
                    case 3
                    % movement by a third -> up or down by a third
                        if this_scale_ind + 2 > 7
                            up_shift = (scale_adjust(this_scale_ind - 5) + 12) - scale_adjust(this_scale_ind);
                        else
                            up_shift = scale_adjust(this_scale_ind + 2) - scale_adjust(this_scale_ind);
                        end
                        if this_scale_ind - 2 <= 0
                            down_shift = scale_adjust(this_scale_ind) - (scale_adjust(this_scale_ind + 5) - 12);
                        else
                            down_shift = scale_adjust(this_scale_ind) - scale_adjust(this_scale_ind - 2);
                        end
                    case 2
                    % movement by the scale -> up or down a scale interval 
                        if this_scale_ind + 1 > 7
                            up_shift = (scale_adjust(this_scale_ind - 6) + 12) - scale_adjust(this_scale_ind);
                        else
                            up_shift = scale_adjust(this_scale_ind + 1) - scale_adjust(this_scale_ind);
                        end
                        if this_scale_ind - 1 <= 0
                            down_shift = scale_adjust(this_scale_ind) - (scale_adjust(this_scale_ind + 6) - 12);
                        else
                            down_shift = scale_adjust(this_scale_ind) - scale_adjust(this_scale_ind - 1);
                        end
                    case 1
                    % movement by a semtione -> up or down by a semitone
                        up_shift = 1;
                        down_shift = 1;
                    case 0
                    % no movement case
                        up_shift = 0;
                        down_shift = 0;
                    otherwise
                    % if note movement is -1
                        up_shift = nan;
                        down_shift = nan;
                        
                end

                % disp(["i = ", i])
                % disp(["delay scaling = ", delay_scaling(i)])
                % disp(["up shift", up_shift])
                % disp(["up prob", up_degree])
                % disp(["down shift", down_shift])
                % disp(["down prob", down_degree])
                % disp(["this movement = ", this_movement])

                if isnan(up_degree)
                % this is just to catch any crazy edge cases
                up_degree = 0.5;
                down_degree = 0.5;
                end
              
                if this_movement ~= -1
                    if (13 + up_shift) > length(doubleside_dist)
                        up_shift = 0;
                    end
                    if (13 - down_shift) < 1
                        down_shift = 0;
                    end

                    % this also encompasses stationary transitions, though the side effect is that these are doubly weighted
                    doubleside_dist(13 + up_shift) = doubleside_dist(13 + up_shift) + (delay_scaling(i) * up_degree);
                    doubleside_dist(13 - down_shift) = doubleside_dist(13 - down_shift) + (delay_scaling(i) * down_degree);
                end

                % % plotting for debugging and algorithm progress representation
                % figure
                % plot(doubleside_dist)
                % hold on
                % plot(last_plot_dist)
                % 
                % pause
                % close all;

                last_plot_dist = doubleside_dist;          
            end

            % now carry out extrapolation to the other octaves
            infer_centre = obj.last_note;

            

            full_bayes_infer = zeros(1,length(obj.note_names));

            upper_lim = infer_centre + 12;
            lower_lim = infer_centre - 12;
            oct = 0;

            for i = infer_centre:length(obj.note_names)
                oct_scaling = floor((i - infer_centre)/13);
                full_bayes_infer(i) = doubleside_dist(mod(12 + (i - infer_centre), 25) + 1) * (1/(oct_falloff ^ oct_scaling));
            end
            for i = infer_centre:-1:1
                oct_scaling = floor((infer_centre - i)/13);
                full_bayes_infer(i) = doubleside_dist(mod(12 - (infer_centre - i), 25) + 1) * (1/(oct_falloff ^ oct_scaling));
            end

            % normalise by mean
            full_bayes_infer = full_bayes_infer ./ mean(full_bayes_infer);
            out_dist = full_bayes_infer;

            % figure
            % plot(full_bayes_infer);
            % xticks(1:length(obj.note_names))
            % xticklabels(obj.note_names(1:end));
            % ylabel("Weighting Assigned to Note")
            % xlabel("Note Names")
            % title("Example of the Output from Bayesian Inference Process")
            % subtitle("Last Note = D3")
            % print(gcf, '-depsc', 'bayes_infer.eps');


            % hold on
            % plot(mean_bayes_infer)

        end
            
    end
end

