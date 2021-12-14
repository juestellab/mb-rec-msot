function [data_crop] = crop_first_n_signals(data_raw,n_signals_to_crop)
    data_crop = data_raw(n_signals_to_crop+1 : end, :, :, :, :);
end
