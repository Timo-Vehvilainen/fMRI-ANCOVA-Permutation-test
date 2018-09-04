clear all
close all

%% Preparations

%load the research data from the repository
memMaps_c1 = load(['/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup' ...
                  '/Resultsinoutgroup/Results/preperation_ISC/' ...
                  'NewAllsubjects_iscMovie1/memMaps']);
memMaps_c1.memMaps.cormatMap.whole.band0.Session1.cor.Filename = ...
    strrep(memMaps_c1.memMaps.cormatMap.whole.band0.Session1.cor.Filename,...
    'Results','Resultsinoutgroup/Results');
memMaps_c2 = load(['/m/nbe/scratch/braindata/afadilm1/Ingroup-outgroup/'
                    'Resultsinoutgroup/Results/preperation_ISC/'
                    'NewAllsubjects_iscMovie2/memMaps']);
memMaps_c2.memMaps.cormatMap.whole.band0.Session1.cor.Filename = ...
    strrep(memMaps_c2.memMaps.cormatMap.whole.band0.Session1.cor.Filename, ...
    'Results','Resultsinoutgroup/Results');

%Loading a mask for selecting the brain out of the nifti-image
mask=load_nii('MNI152_T1_2mm_brain_mask.nii');


%% The main loop



%looping through all the voxels in the nifti-image
parfor x = 1:91
    
    %Initiating the current slice of voxels.
    %The results are saved after every slice in case of a malfunction.
    X_slice_p = zeros(109, 91);
    X_slice_t = zeros(109, 91);
    
    for y = 1:109
        for z = 1:91
            
            %only process this voxel if it is inside the brain, otherwise
            %skip it and go to the next voxel
            if(mask.img(x,y,z)==0)             
                continue; 
            end
            
            %Call the function for calculating the p-value, and store it
            [X_slice_p(y, z) X_slice_t(y, z)] = ...
                parfunction(memMaps_c1, memMaps_c2, x, y, z);
        end
    end
    %Save this 2d-slice of finished voxels in a file
    parsave(x, X_slice_p);
end

%% Combining the slices into a full 3D voxel matrix
% 
 result_g1vsg2_p = zeros(91, 109, 91);
 result_g1vsg2_t = zeros(91, 109, 91);
 for i = 1:91
     slice_t = load([pwd '/tstats/X_slice' num2str(i) '.mat']);
     slice_p = load([pwd '/slices/X_slice' num2str(i) '.mat']);
     result_g1vsg2_p(i, :, :) = slice_p.X_slice;
     result_g1vsg2_t(i, :, :) = slice_t.X_slice;
 end
 save('result_g1vsg2_tstat.mat', 'result_g1vsg2_t');
 save('result_g1vsg2_pvals.mat', 'result_g1vsg2_p');
 
 %% Making p-values into q-values via FDR-correction
 
%Forming the p-values into a vector for the mafdr()-function
pvector = reshape(result_g1vsg2_p, [1 902629]);

%Calculating the q-values
q=mafdr(pvector,'BHFDR','true');

%separating the significant indices and making a mask of them
q_unsignificant_indices = find(q >= 0.05);
q_significant_indices = find(q < 0.05);
q_mask = q;
q_mask(q_unsignificant_indices) = 0;
q_mask(q_significant_indices) = 1;
q_mask = reshape(q_mask, [91, 109, 91]);

%We use the mask to filter out FDR-corrected t-scores
result_g1vsg2_t_qmasked = result_g1vsg2_t;
result_g1vsg2_t_qmasked(q_mask == 0) = 0;

%Save the resulting data in a nifti-image
result_g1vsg2_t_qmasked_nii = make_nii(result_g1vsg2_t_qmasked);
save_nii(result_g1vsg2_t_qmasked_nii, 'result_g1vsg2_t_qmasked.nii');
result_g1vsg2_t_qmasked_fixedOriginator = ...
    bramila_fixOriginator('result_g1vsg2_t_qmasked.nii');
save_nii(result_g1vsg2_t_qmasked_fixedOriginator, ...
    'result_g1vsg2_t_qmasked_fixedOriginator.nii');
