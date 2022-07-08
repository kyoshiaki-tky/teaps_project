function [dist_conv,Teaps_output_all]  = dist_conv_judge(Teaps_output,alpha)
%DIST_CONV_JUDGE 

% retrieve iteratioln number of whole teaps 
iterN = size(Teaps_output,2);
freeXN = size(Teaps_output{1}.finalX,1);

Teaps_output_all =struct;
Teaps_output_all.finalX = zeros(size(Teaps_output{1}.finalX,1),...
    size(Teaps_output{1}.finalX,2)*iterN);

Teaps_output_prev = zeros(size(Teaps_output{1}.finalX,1),...
    size(Teaps_output{1}.finalX,2)*(iterN-1));

for k = 1:iterN-1
    start_idx = 1+(k-1)*size(Teaps_output{k}.finalX,2);
    end_idx = k*size(Teaps_output{k}.finalX,2);
    Teaps_output_prev(:,start_idx:end_idx) = Teaps_output{k}.finalX;
end

Teaps_output_all.finalX(:,:) = [Teaps_output_prev, Teaps_output{iterN}.finalX];

% wilcoxon rank sum test against each row
p_vector = zeros(1,freeXN);
for k = 1:freeXN
    p_vector(k) = ranksum(Teaps_output_prev(k,:),Teaps_output_all.finalX(k,:));
end

dist_conv = (sum(p_vector<alpha)==0);

end

