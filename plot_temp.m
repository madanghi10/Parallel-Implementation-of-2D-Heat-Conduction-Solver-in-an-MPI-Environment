clear
tid = [0,1052,2104,2631];

ranks = 0; 

px = 1; 
py = 1; 

for i = 1:length(tid)
    all_data = [];
    proc_dims = zeros(length(ranks), 4); 
    
    for r = 1:length(ranks)
        rank = ranks(r);
        filename = sprintf('T_x_y_%06d_%04d_%d*%d.dat', tid(i), rank,px,py);
        data = dlmread(filename);
        
        x_vals = unique(data(:,1));
        y_vals = unique(data(:,2));
        nx = length(x_vals);
        ny = length(y_vals);
        
        proc_dims(r,:) = [rank, nx, ny, size(all_data,1)+1];
        all_data = [all_data; data];
    end
    
    x_global = unique(all_data(:,1));
    y_global = unique(all_data(:,2));
    nx_global = length(x_global);
    ny_global = length(y_global);
    T_global = zeros(nx_global, ny_global);
    
    for r = 1:length(ranks)
        rank = proc_dims(r,1);
        nx = proc_dims(r,2);
        ny = proc_dims(r,3);
        start_idx = proc_dims(r,4);
        
        data = all_data(start_idx:start_idx+nx*ny-1,:);
        rank_x = mod(rank, px);
        rank_y = floor(rank / px);
        x_local = unique(data(:,1));
        y_local = unique(data(:,2));
        
        [~, x_start] = min(abs(x_global - x_local(1)));
        [~, y_start] = min(abs(y_global - y_local(1)));
        
        T_local = reshape(data(:,3), [ny, nx])'; 
        T_global(x_start:x_start+nx-1, y_start:y_start+ny-1) = T_local;
    end
    figure;
    [X,Y] = meshgrid(x_global, y_global);
    contourf(X, Y, T_global', 'LineColor', 'none');
    xlabel('x'); ylabel('y');
    title(sprintf('At timestep = %06d using p=%d*%d', tid(i),px,py));
    colormap('jet')
    colorbar;
    saveas(gcf, sprintf('cont_T_%06d_%d_%d.png', tid(i),px,py));
    
end