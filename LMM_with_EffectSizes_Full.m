function results = LMM_with_EffectSizes_Full(tb, formula, x_list, y_list)
    all_needed_vars = unique([x_list(:)', y_list(:)', extractVarsFromFormula(formula)]);
    for v = 1:length(all_needed_vars)
        v_name = all_needed_vars{v};
        if iscell(tb.(v_name))
            if ischar(tb.(v_name){1}) || isstring(tb.(v_name){1})
                tb.(v_name) = str2double(tb.(v_name));
            else
                tb.(v_name) = cell2mat(tb.(v_name));
            end
        end
    end
    % -----------------------------------------------

    total_rows = length(x_list) * length(y_list);
    results_cell = cell(total_rows, 1);
    
    parfor idx = 1:total_rows
        [j, k] = ind2sub([length(y_list), length(x_list)], idx);
        y_name = y_list{j};
        x_name = x_list{k};
        curr_tb = tb(:, unique([y_name, x_name, extractVarsFromFormula(formula)]));
        curr_tb.(y_name) = (curr_tb.(y_name) - nanmean(curr_tb.(y_name))) / nanstd(curr_tb.(y_name));
        is_binary_x = false;
        x_data = curr_tb.(x_name);
        if iscell(x_data)
            valid_idx = ~cellfun(@(x) any(isnan(x)), x_data);
        else
            valid_idx = ~isnan(x_data);
        end
        unique_vals = unique(x_data(valid_idx));
        if length(unique_vals) == 2 && all(ismember(unique_vals, [0, 1]))
            is_binary_x = true;
        else
            curr_tb.(x_name) = (x_data - nanmean(x_data)) / nanstd(x_data);
        end
        
        % Covariates
        covars = extractVarsFromFormula(formula);
        for v = 1:length(covars)
            v_name = covars{v};
            v_data = curr_tb.(v_name);
            if isnumeric(v_data)
                u_vals = unique(v_data(~isnan(v_data)));
                if ~(length(u_vals) == 2 && all(ismember(u_vals, [0, 1])))
                    curr_tb.(v_name) = (v_data - nanmean(v_data)) / nanstd(v_data);
                end
            end
        end

        formula_input = [y_name, ' ~ ', x_name, ' + ', formula];
        
        try
            lmm = fitlme(curr_tb, formula_input);
            fe_table = lmm.Coefficients;
            target_row = fe_table(strcmp(fe_table.Name, x_name), :);

            if ~isempty(target_row)
                t = target_row.tStat;
                p = target_row.pValue;
                DF = target_row.DF;
                std_beta = target_row.Estimate; 
                ci1 = target_row.Lower;
                ci2 = target_row.Upper;

                % Effect Sizes
                partial_R2 = t^2 / (t^2 + DF);

                results_cell{idx} = {x_name, y_name, t, p, DF, std_beta, ci1, ci2, is_binary_x, ...
                                     partial_R2};
            end
            disp(['Linear Mixed Model: ', y_list{j}, ' ~ ', x_list{k}, ' ']);
            
        catch ME
            results_cell{idx} = {x_name, y_name, NaN, NaN, NaN, ...
                                 NaN, NaN, NaN, is_binary_x, NaN};
        end
    end

    results = cell2table(vertcat(results_cell{:}), ...
        'VariableNames', {'xname', 'yname', 't', 'p', 'DF', 'std_beta', 'ci1', 'ci2', 'is_binary', ...
                          'Partial_R2'});
end

function vars = extractVarsFromFormula(f)
    tokens = regexp(f, '\w+', 'match');
    keywords = {'1', 'Site', 'Subject', 'ID', '2'}; 
    vars = tokens(~ismember(tokens, keywords));
end
