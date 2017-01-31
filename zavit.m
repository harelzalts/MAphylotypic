function x = zavit(y);

%%%x is the gene expression order indices

%%%This function enables to sort temporal gene expression data from early
%%%to late. This is done in a few basic stages: First, we standadize the 
%%%data so that for each row, which represents a temporal expression vector
%%%of a specific gene, the mean is 0 with a standard deviation of 1. Then,
%%%we preform PC analysis on this data, and look only at PC1 and PC2. 
%%%The PC1 and PC2 scatter plot creates a circle. We proceed to calculate
%%%the inverse tangent for each point (which represents a single gene)
%%%in the scatter plot, then sort these values to obtain a temporal order.


%%% 1. Standardize the data (using the function zscore)
a = zscore(y,0,2);

%%% 2. Principal Component analysis of the standardized data
[pc, scores] = pca(a);

%%% 3. Calculate the inverse tangent of each gene, then sort the genes
%%% according to the obtained angle 
ZAVIT = mod(atan2d(scores(:,1),scores(:,2)), 0);
[~, x] = sort(ZAVIT);

corr_vector = [1 repmat(0,1,size(a,2)-1)];

%%%Choose the first gene
correlation_with_start = corrcoef([corr_vector; a(x,:) ]');
[m, i_m] = max(correlation_with_start(1,2:end));
x = [x(i_m:end)', x(1:i_m-1)'];

[m, m_i]= max(a(x,:)');

c = diag(corrcoef(m_i, 1:length(m_i)),1);

if (c<0)
    x = flipud(x');
end

end