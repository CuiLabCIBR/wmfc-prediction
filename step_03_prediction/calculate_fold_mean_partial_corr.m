function [mean_partial_gw, mean_partial_ww, fold_partial_gw, fold_partial_ww] = calculate_fold_mean_partial_corr(base_folder, repetition_index, num_folds)
mean_partial_gw = NaN;
mean_partial_ww = NaN;
fold_partial_gw = nan(1, num_folds);
fold_partial_ww = nan(1, num_folds);

for fold_index = 0:(num_folds - 1)
    gg_file = fullfile(base_folder, ['Time_' num2str(repetition_index)], 'GGFC', ['Fold_' num2str(fold_index) '_Score.mat']);
    gw_file = fullfile(base_folder, ['Time_' num2str(repetition_index)], 'GWFC', ['Fold_' num2str(fold_index) '_Score.mat']);
    ww_file = fullfile(base_folder, ['Time_' num2str(repetition_index)], 'WWFC', ['Fold_' num2str(fold_index) '_Score.mat']);

    if ~isfile(gg_file) || ~isfile(gw_file) || ~isfile(ww_file)
        warning('Missing fold score file for repetition %d, fold %d. Skipping this repetition.', repetition_index, fold_index);
        return;
    end

    try
        gg = load(gg_file, 'Index', 'Predict_Score', 'Test_Score');
        gw = load(gw_file, 'Index', 'Predict_Score', 'Test_Score');
        ww = load(ww_file, 'Index', 'Predict_Score', 'Test_Score');

        [gg_index, gg_order] = sort(double(gg.Index(:)));
        [gw_index, gw_order] = sort(double(gw.Index(:)));
        [ww_index, ww_order] = sort(double(ww.Index(:)));

        if ~isequal(gg_index, gw_index) || ~isequal(gg_index, ww_index)
            warning('Participant indices differ across FC types for repetition %d, fold %d. Skipping this repetition.', repetition_index, fold_index);
            return;
        end

        gg_predict = double(gg.Predict_Score(:));
        gg_test = double(gg.Test_Score(:));
        gw_predict = double(gw.Predict_Score(:));
        gw_test = double(gw.Test_Score(:));
        ww_predict = double(ww.Predict_Score(:));
        ww_test = double(ww.Test_Score(:));

        gg_predict = gg_predict(gg_order);
        gg_test = gg_test(gg_order);
        gw_predict = gw_predict(gw_order);
        gw_test = gw_test(gw_order);
        ww_predict = ww_predict(ww_order);
        ww_test = ww_test(ww_order);

        if ~isequaln(gg_test, gw_test) || ~isequaln(gg_test, ww_test)
            warning('Test scores differ across FC types for repetition %d, fold %d. Skipping this repetition.', repetition_index, fold_index);
            return;
        end

        fold_partial_gw(fold_index + 1) = partialcorr(gw_predict, gw_test, gg_predict);
        fold_partial_ww(fold_index + 1) = partialcorr(ww_predict, ww_test, gg_predict);
    catch ME
        warning('Partial correlation failed for repetition %d, fold %d: %s. Skipping this repetition.', repetition_index, fold_index, ME.message);
        fold_partial_gw(:) = NaN;
        fold_partial_ww(:) = NaN;
        return;
    end
end

if all(isfinite(fold_partial_gw)) && all(isfinite(fold_partial_ww))
    mean_partial_gw = mean(fold_partial_gw);
    mean_partial_ww = mean(fold_partial_ww);
else
    fold_partial_gw(:) = NaN;
    fold_partial_ww(:) = NaN;
end
end
