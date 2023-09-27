function retour = Cizdir9_Newmod(array,idx,idx_array,size_idx)
% Change the order of row of ARRAY, starting by the row number IDX

indexes = idx_array(idx:(idx+size_idx-1));
retour = array(indexes,:);