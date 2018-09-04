function parsave(x, X_slice_p, X_slice_t )
%PARSAVE Function for saving the finished slices into .mat-files in a
%parallizabe for-loop
    save([pwd '/slices/X_slice' num2str(x) '.mat'], 'X_slice_p');
    save([pwd '/tstats/X_slice' num2str(x) '.mat'], 'X_slice_t');
end

