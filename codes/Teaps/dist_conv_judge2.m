function dist_conv  = dist_conv_judge2(Teaps_output,beta)
%DIST_CONV_JUDGE2 

N_iter=size(Teaps_output,2);
for k = 1:(N_iter-1)
    if k==1
        max_vec = max(Teaps_output{1,k}.finalX,[],2);
        min_vec = min(Teaps_output{1,k}.finalX,[],2);
    else
        max_vec_tmp = max(Teaps_output{1,k}.finalX,[],2);
        min_vec_tmp = min(Teaps_output{1,k}.finalX,[],2);
        max_vec = max([max_vec_tmp, max_vec],[],2);
        min_vec = min([min_vec_tmp, min_vec],[],2);
        clear max_vec_tmp min_vec_tmp
    end
end

n_outside =0;
Teaps_output{1,N_iter}.finalX;
N_sample = size(Teaps_output{1,N_iter}.finalX,2);
for k = 1:N_sample
    % count max and min violation
    if sum(Teaps_output{1,N_iter}.finalX(:,k)>max_vec)>0
        n_outside=n_outside+1;
    elseif sum(Teaps_output{1,N_iter}.finalX(:,k)<min_vec)>0
        n_outside=n_outside+1;
    else
    end
end

inside_ratio = 1-n_outside/N_sample;
dist_conv = (inside_ratio >= beta);

end

