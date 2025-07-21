clear
tid = 10;

ranks = 0:15; 

px = 4; 
py = 4; 
serial_data=dlmread("T_x_y_000010_0000_1*1.dat");
n_s=sqrt(size(serial_data,1));
T_serial=reshape(serial_data(:,3),[n_s,n_s]);
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
end
compute_error = @(A, B) deal(sqrt(mean((A - B).^2, 'all')));

[L2_p2_2] = compute_error(T_serial, T_global); % Serial vs p=2*2
fprintf('Error (serial vs p=2): L2 = %350e\n', L2_p2_2);