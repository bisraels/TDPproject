%Header Goes Here
%--------------------------------------------------------------------------
% Find local maxima and values of peaks
%--------------------------------------------------------------------------
function [Num_loc_max, mtx_trns] = locMax_finder(tdpGrid, tdpGrid_gauss_normalized)

L = length(tdpGrid);
Z = tdpGrid_gauss_normalized;   % Z-axis of the surface plot

Num_loc_max = 0;

%mtx_trns = zeros(1,3);

for i = 2:(L-1)
    for j = 2:(L-1)
        
        cell = Z(i,j);          % Z is the matrix that we are plotting in the surface plot.
        
        % Evaluate neighboring matrix cells
        above = Z(i, j + 1);
        below = Z(i, j - 1);
        left = Z(i - 1, j);
        right = Z(i + 1,j);
        
        if above < cell && below < cell && left < cell && right < cell
            Num_loc_max = Num_loc_max + 1;
           % disp('Local max found at (x ,y, z) =');
           % disp([(i)/201,(j)/201, Z(i,j)]);     
            
            if exist('mtx_trns')
                mtx_trns = vertcat(mtx_trns, [i/gridRes,j/gridRes, Z(i,j)]);
            else
                mtx_trns = [(i)/gridRes,(j)/gridRes, Z(i,j)];
            end
        end
    end
end

Num_loc_max;
mtx_trns;

% Plot a red dot on each local maxima accounted for:
plot3(mtx_trns(:,1),mtx_trns(:,2),mtx_trns(:,3),'r.','MarkerSize',14)

end