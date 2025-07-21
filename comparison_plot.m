clear;
tid = 2631; %change the value
px = [1, 2, 2, 4];
py = [1, 2, 4, 4];

figure;  
hold on; 

legend_entries = {}; 

for i = 1:length(tid)  
    for j = 1:length(px)  
        ranks = 0:((px(j) * py(j)) - 1);
        all_data = [];
        proc_dims = zeros(length(ranks), 4);

        
        for r = 1:length(ranks)
            rank = ranks(r);

            
            filename = sprintf('T_x_y_%06d_%04d_%d*%d.dat', tid(i), rank, px(j), py(j));
            data = dlmread(filename);
            
       
            x_local = unique(data(:, 1));
            y_local = unique(data(:, 2));
            
            
            nx = length(x_local);
            ny = length(y_local);
            
            proc_dims(r, :) = [rank, nx, ny, size(all_data, 1) + 1];
            
            all_data = [all_data; data];
        end
        
        x_global = unique(all_data(:, 1));
        y_global = unique(all_data(:, 2));
        nx_global = length(x_global);
        ny_global = length(y_global);
        
        T_global = zeros(nx_global, ny_global);

        for r = 1:length(ranks)
            rank = proc_dims(r, 1);
            nx = proc_dims(r, 2);
            ny = proc_dims(r, 3);
            start_idx = proc_dims(r, 4);
            
            data = all_data(start_idx:start_idx + nx * ny - 1, :);

            x_vals = unique(data(:, 1));
            y_vals = unique(data(:, 2));
            
            [~, x_start] = min(abs(x_global - x_vals(1)));
            [~, y_start] = min(abs(y_global - y_vals(1)));

            T_local = reshape(data(:, 3), [ny, nx])';
            
            T_global(x_start:x_start + nx - 1, y_start:y_start + ny - 1) = T_local;
        end
        
        [~, mid_idx] = min(abs(y_global - 0.5));
        plot(x_global, T_global(:, mid_idx), '-', 'LineWidth', 2);

        legend_entries{end + 1} = sprintf('tid=%d, %dx%d', tid(i), px(j), py(j));
    end
end

xlabel('x');
ylabel('T');
title(sprintf('Comparison of Profile at y=0.5, timestep = %06d ', tid(i)));
legend(legend_entries);

saveas(gcf, 'Comparison_Profile_tid_2631.png');%change tid value while saving
