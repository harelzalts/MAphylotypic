
clear all
close all

try
    values_tpm_log;
catch    
    load Raw_data.mat;
    load Gene_names.mat;
    load samples_and_strain_names.mat
    load Normalized_data.mat;
    load corr_matrix.mat;    
    load MA_bio_replicate_raw.mat;
    load MA_bioreplicate_norm;
    load distance_RW_data;
    load distance_RW_replicates_data;
    load clean_distance.mat;
    load Gene_names_RW;    
    load TH_lists.mat;
    load special_gene_lists.mat;
    load sig_paths2;
    load TF_matrix;
    load GA_matrix.mat;
    load GO_tables.mat
    load cmap_color_blind2.mat
    load tpm_50.mat;
    load whole_tc_tpm;
    load running_window_expression_pattern;
    
end


preperation_normalization                           = 1;

choose_random_quad                                  = 0;

General_stuff                                       = 1;

supplementary_figure_1                              = 0;

Calculate_distance                                  = 1;

running_window_expression_pattern                   = 0;

replicates_analysis                                 = 0;

distance_cleanup_and_running_window                 = 0;

GO_analysis                                         = 0;

functional_groups_analysis                          = 0;

figure_4_analysis                                   = 0;

control_analysis                                    = 0;

GO_enrichment_for_gene_subgroups                    = 0;

RW_to_stage_sorting                                 = 0;


% Parameters:

log_norm                    = 1;
manual_move_genes_in_ZAVIT  = 1;
permutation_ind             = -5410;

number_of_experiment_points = 7;
number_of_repeats           = 20;

window_size                 = 200;
step_size                   = 100;


if (preperation_normalization)
         
    %Additional UMI normalization
    values = real(-(4^5)*log(1-values/(4^5)));
    
    %TPM normalization (actually, TP per 50,000)
    values_tpm = zeros(size(values));
    for i=1:size(values,2)
        values_tpm(:,i) = values(:,i)*50000/(sum(values(:,i),1));
    end
    
    %log normalization (log2)
    if (log_norm)
        values_tpm_log = log2(1+values_tpm);
    else
        values_tpm_log = values_tpm;
    end
    
    %remove lincs
    values_tpm_log = values_tpm_log(1:20517,:);
       
    %define strain names
    strains = {'MA 550' 'MA 570' 'MA 502' 'MA 555' 'N2'...
        'MA 575' 'MA 530' 'MA 517' 'MA 529' 'MA 598' 'MA 593'...
        'MA 573' 'MA 564' 'MA 516' 'MA 552' 'MA 522'...
        'MA 536' 'MA 508' 'MA 566' 'MA 597'};
    
    %Create work matrices:
    %MALib_TP = Expression divided by time points, each time point in a
    %different matrix
    %all_time_points = Expression divided by time points in one matrix
    %mean_expression_matrix = the mean expression of 20 strains in each
    %time point
    %strains_in_3D = Expression divided by strains, each strain in a
    %different matrix
    
    MALib_TP = [];
    sample_names = {};
    for i = 1:number_of_experiment_points
        MALib_TP(:,:,i) = values_tpm_log(:, i:7:size(values_tpm_log,2));
        sample_names(:,i) = Sample_names(1,i:7:size(Sample_names,2));
    end
    
    all_time_points = [MALib_TP(:,:,1), MALib_TP(:,:,2),MALib_TP(:,:,3),...
        MALib_TP(:,:,4),MALib_TP(:,:,5),MALib_TP(:,:,6),MALib_TP(:,:,7)];
    samples = vertcat(sample_names(:,1),sample_names(:,2),sample_names(:,3),...
        sample_names(:,4),sample_names(:,5),sample_names(:,6),sample_names(:,7));
    
    mean_expression_matrix = []; 
    for i = 1:number_of_experiment_points
        mean_expression_matrix(:,i) = mean(MALib_TP(:,:,i),2);
    end
    
    strains_in_3D = [];
    for i = 1:number_of_repeats
        strains_in_3D(:,:,i) = all_time_points(:,i:20:size(all_time_points,2));
    end
         
    try
        load corr_matrix.mat;
    catch
              
        corr_matrix = [];
        for i = 1:number_of_repeats
            corr_matrix(:,i) = diag(corr(mean_expression_matrix',strains_in_3D(:,:,i)'));
        end
        meanCorr = mean(corr_matrix,2);
               
        
        save corr_matrix.mat corr_matrix meanCorr
    end
           
    save Normalized_data.mat ...
        values_tpm_log mean_expression_matrix...
        all_time_points strains_in_3D MALib_TP;

end

if (choose_random_quad)
    
    d = randperm(20);
    rand_quad = d(1:4);
    values_tpm_rand = [];
    for i = 1:length(rand_quad)
        for l = rand_quad(1,i)
            values_tpm_rand(:,:,i) = all_time_points(:,[l:number_of_repeats:size(values,2)]);
            strains_in_3D_rand(:,:,i) = strains_in_3D(:,:,rand_quad(i));
        end
    end
    
    values_tpm_log_rand = [values_tpm_rand(:,:,1) values_tpm_rand(:,:,2) values_tpm_rand(:,:,3) values_tpm_rand(:,:,4)];
            
    MALib_TP_rand = [];
    for i = 1:number_of_experiment_points
        MALib_TP_rand(:,:,i) = values_tpm_log_rand(:, i:7:size(values_tpm_log_rand,2));
    end
    
    all_time_points_rand = [MALib_TP_rand(:,:,1), MALib_TP_rand(:,:,2),MALib_TP_rand(:,:,3),...
        MALib_TP_rand(:,:,4),MALib_TP_rand(:,:,5),MALib_TP_rand(:,:,6),MALib_TP_rand(:,:,7)];
    
    mean_expression_matrix_rand = []; %mean expression matrix
    for i = 1:number_of_experiment_points
        mean_expression_matrix_rand(:,i) = mean(MALib_TP_rand(:,:,i),2);
    end
    
  
    values_tpm_log = [];all_time_points = [];all_time_points = [];MALib_TP = [];strains_in_3D = [];
    mean_expression_matrix = [];corr_matrix = [];
    
    values_tpm_log = values_tpm_log_rand;
    
    all_time_points = all_time_points_rand;
    
    MALib_TP = MALib_TP_rand;
    
    strains_in_3D = strains_in_3D_rand;
    
    mean_expression_matrix = mean_expression_matrix_rand;
    

       
end

if (General_stuff)%(Figure 1B,C and 2B)
    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
       
    strains_in_3D_2 = strains_in_3D;
    strains_in_3D_2(:,:,1) = strains_in_3D(:,:,5);
    strains_in_3D_2(:,:,5) = strains_in_3D(:,:,1);
    
    strains_a = {'N2' 'MA 570' 'MA 502' 'MA 555' 'MA 550' ...
        'MA 575' 'MA 530' 'MA 517' 'MA 529' 'MA 598' 'MA 593'...
        'MA 573' 'MA 564' 'MA 516' 'MA 552' 'MA 522'...
        'MA 536' 'MA 508' 'MA 566' 'MA 597'};
    
    all_time_points_2 = all_time_points;
    all_time_points_2(1:20:140) = all_time_points(5:20:140);
    all_time_points_2(5:20:140) = all_time_points(1:20:140);
    

    figure;imagesc(corr(all_time_points_2(selected_genes,:)),[0 1]);colormap(cmap_color_blind2);colorbar
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
        
    
    %%plots for single genes

    cmap = colormap(jet(size(strains_in_3D,3)));
    
        for k =  find(ismember(Genes_no_lincs(selected_genes) ,'T26E3.3'));
            j = selected_genes(k);
            
            figure;
            for i = 1:size(strains_in_3D_2,3)
%                 plot(strains_in_3D_2(j,:,i),'LineWidth',3,'color',cmap(i,:))
                plot(smooth(strains_in_3D_2(j,:,i)),'LineWidth',3,'color',cmap(i,:))
                hold on
            end
            title(CE_official_gene_name(j),'FontSize', 20, 'FontWeight', 'bold')
%             xlabel('Time', 'FontSize', 16, 'FontWeight', 'bold')
%             ylabel('Expression Level (Log2TPM)', 'FontSize', 16, 'FontWeight', 'bold')
%             legend(strains_a)
            set(gca,'FontSize',16,'FontWeight','bold');
%             set(gca,'xtick',[],'ytick',[]);
%             xlim([1 7]);ylim([-0.1 3.5])
        end
    
        
    %%This is for figure2b - heat map of gene indices
   
    load gene_pos_ind_replicates;load gene_pos_ind_in_means_matrix
    
    gene_pos_ind_ref = (1:length(selected_genes))';
    blat = [gene_order_reference gene_order_each_strain];

    [dd,ddd] = sort(mean_distance_temporal_clean,'ascend');
    sort_genes_by_distnace = gene_names_order_distance_matrix(ddd);
    
    gene_of_interest1 = {'F44C4.5'};
    gene_of_interest2 = {'C25A1.13'};

    gene_of_interest_in_each_strain1 = [];gene_of_interest_in_each_strain2 = [];
    for i = 1:21;
        gene_of_interest_in_each_strain1(i) = find(ismember(...
            blat(:,i),gene_of_interest1));
        gene_of_interest_in_each_strain2(i) = find(ismember(...
            blat(:,i),gene_of_interest2));
    end
    
    figure;
    subplot(1,2,2);
    h = imagesc([gene_pos_ind_ref gene_pos_ind_each_strain]);colormap(paruly);colorbar;
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    hold on
    plot(gene_of_interest_in_each_strain1,'*black')
    plot(gene_of_interest_in_each_strain2,'*magenta')


    blat2 = [Genes_sorted gene_order_each_replicate];
    gene_of_interest_in_each_strain_reps1 = [];gene_of_interest_in_each_strain_reps2 = [];
    for i = 1:size(blat2,2)
        gene_of_interest_in_each_strain_reps1(i) = find(ismember(...
            blat2(:,i),gene_of_interest1));
        gene_of_interest_in_each_strain_reps2(i) = find(ismember(...
            blat2(:,i),gene_of_interest2));
    end
    
    subplot(1,2,1);
    imagesc([gene_pos_ind_ref gene_pos_ind]);colormap(paruly)
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    hold on
    plot(gene_of_interest_in_each_strain_reps1,'*black')
    plot(gene_of_interest_in_each_strain_reps2,'*magenta')
      
end

if (supplementary_figure_1)
    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
    
   % Mean Corrrelation of each time point transcriptome
    
    j = 1:20:140;
    if (choose_random_quad)
        j = 1:4:28;
    end
    corr_TP = {i};mean_corr_TP = [];corr_TP_2 = [];
    for i = 1:7
        if (choose_random_quad == 0)
            corr_TP{i} = corr(all_time_points(selected_genes,j(i):j(i)+19));
        elseif (choose_random_quad == 1)
            corr_TP{i} = corr(all_time_points(selected_genes,j(i):j(i)+3));
        end
        corr_TP_2(:,i) = corr_TP{i}(find(triu(corr_TP{i},1)));
        mean_corr_TP(i) = mean(corr_TP{i}(find(triu(corr_TP{i},1))));
        std_corr_TP(i) = std(corr_TP{i}(find(triu(corr_TP{i},1))));
    end   
    
    j = 1:4:28;
    corr_TP_br = {};mean_corr_TP_br = [];corr_TP_br_2 = [];
    for i = 1:7
        corr_TP_br{i} = corr(biorepeat_expression_tpm_log(selected_genes,j(i):j(i)+3));
        
        corr_TP_br_2(:,i) = corr_TP_br{i}(find(triu(corr_TP_br{i},1)));
        
        mean_corr_TP_br(i) = mean(corr_TP_br{i}(find(triu(corr_TP_br{i},1))));
        std_corr_TP_br(i) = std(corr_TP_br{i}(find(triu(corr_TP_br{i},1))));
    end
            
    if (choose_random_quad)
        strains(rand_quad)
    end
    
    corr_TP_br_2(7:size(corr_TP_2,1),:) = nan;
      
    corr_TP_3 = [corr_TP_2(:,1) corr_TP_br_2(:,1) corr_TP_2(:,2) corr_TP_br_2(:,2)...
        corr_TP_2(:,3) corr_TP_br_2(:,3) corr_TP_2(:,4) corr_TP_br_2(:,4)...
        corr_TP_2(:,5) corr_TP_br_2(:,5) corr_TP_2(:,6) corr_TP_br_2(:,6)...
        corr_TP_2(:,7) corr_TP_br_2(:,7)];
    
    figure;
    h = boxplot(corr_TP_3);ylim([0.5 1])
    set(h([ 7],:),'Visible','off');set(h(5,:),'LineWidth',1)
    set(gca,'ytick',[],'yticklabel',[],'xtick',[],'xticklabel',[]);
    h2 = findobj(gca,'Tag','Box');

    nn = [1:2:13 2:2:14];
    for j=1:7
        patch(get(h2(nn(j)),'XData'),get(h2(nn(j)),'YData'),[0.8 0.6 0.2],'FaceAlpha',.5);
    end;
    for j=8:14
        patch(get(h2(nn(j)),'XData'),get(h2(nn(j)),'YData'),[0.1 0.2 0.6],'FaceAlpha',.5);
    end;
    
end

if (Calculate_distance)
        
    %In this section, we calculate the distance between the expression
    %timing of temporally arranged genes in the 20 MA strains 
    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    
    %Select dynamic genes:
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
        
    %Order the genes using ZAVIT:
    
    %1. Determine the first gene
    
    Z_expression_matrix = zscore(values_tpm_log(selected_genes,:),0,2);
    [pcg, zscoresg] = pca(Z_expression_matrix);
    
    zavit = mod(atan2d(zscoresg(:,1),zscoresg(:,2)), 0);
    [~, gene_order_of_expression_matrix] = sort(zavit);
           
    z_all_time_points = zscore(all_time_points(selected_genes,:),0,2);
    
    figure;imagesc(z_all_time_points(gene_order_of_expression_matrix,:));colormap(cmap_color_blind2);colorbar
    set(gca, 'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    
    figure;imagesc(Z_expression_matrix(gene_order_of_expression_matrix,:));colormap(cmap_color_blind2);colorbar;
    set(gca, 'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
   
    pcg = [];zscoresg = [];pcvarsg = [];zavit = [];gene_order = [];
    each_strain_fitted = [];zavit = [];
    

    %2. ZAVIT for each strain seperately
    
%     try
%         load standardized_data_strains_with_order_by_zavit
%     catch
        for i = 1:size(strains_in_3D,3);
            
            each_strain_fitted(:,:,i) = zscore(strains_in_3D(selected_genes,:,i),0,2);
            
            [pcg(:,:,i), zscoresg(:,:,i)] = pca(each_strain_fitted(:,:,i));
            zavit(:,i) = mod(atan2d(zscoresg(:,1,i),zscoresg(:,2,i)), 0) ;
            [~, gene_order(:,i)] = sort(zavit(:,i));
            
        end
        
        %This produces a messy picture. To fix this, we use the first gene
        %found above
               
        each_strain_fitted_and_sorted = [];gene_reorder = []; 
        for i = 1:size(each_strain_fitted,3)
                      
            first_gene(i) = find(gene_order_of_expression_matrix(1) == gene_order(:,i));

            gene_reorder(:,i) = circshift(gene_order(:,i),(size(gene_order,1)-(first_gene(i)-1)));
            
            [~,m_i] = max(each_strain_fitted(gene_reorder(:,i),:,i));
            c  =corrcoef(m_i, 1:length(m_i)); c = c(1,2);
            if (c<0.3)
                gene_reorder(:,i) = flipud(gene_reorder(:,i));
            end
            
            
        end
%     end
    
    if (manual_move_genes_in_ZAVIT)        
       gene_reorder = circshift(gene_reorder,permutation_ind); 
    end
    
    for i = 1:size(strains_in_3D,3);
        each_strain_fitted_and_sorted(:,:,i) = each_strain_fitted(gene_reorder(:,i),:,i);
    end
    
    figure;
    for i  = 1:20
    subplot(4,5,i)
        imagesc(each_strain_fitted_and_sorted(:,:,i));colormap(redblue)
    end
    
    
    
    save standardized_data_strains_with_order_by_zavit each_strain_fitted...
        gene_reorder each_strain_fitted_and_sorted

    
    %3. ZAVIT for the means matrix
    
    Z_mean_expression_matrix = zscore(mean_expression_matrix(selected_genes,:),0,2);
    [pcg, zscoresg] = pca(Z_mean_expression_matrix);
    
    zavit = mod(atan2d(zscoresg(:,1),zscoresg(:,2)), 0);
    [~, gene_order_of_mean_expression_matrix] = sort(zavit);
    
    figure;
    scatter(zscoresg(:,1),zscoresg(:,2),10,zavit,'fill');colormap(paruly);colorbar(); 
    
    correlation_with_start = corrcoef([ 1 0 0 0 0 0 0;...
        Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:) ]');
    
    [m, i_m] = max(correlation_with_start(1,2:end));
    
    gene_order_of_mean_expression_matrix = [gene_order_of_mean_expression_matrix(i_m:end)', gene_order_of_mean_expression_matrix(1:i_m-1)'];
    
    [m, m_i]=max(Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:)');
    
    c = diag(corrcoef(m_i, 1:length(m_i)),1);
    
    if (c<0)
        gene_order_of_mean_expression_matrix = flipud(gene_order_of_mean_expression_matrix');
    end
    
    if (manual_move_genes_in_ZAVIT)        
       gene_order_of_mean_expression_matrix = circshift(gene_order_of_mean_expression_matrix',permutation_ind); 
    end
    
    figure;imagesc(Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:));colormap(redblue);colorbar;
    set(gca, 'xtick',[],'xticklabel',[],'ytick',[1000 2000 3000 4000 5000],...
        'yticklabel',[],'TickDir','out','TickLength',[0.03 0]);
 

    %     %%%now that I have the gene order in each strain, I create another
    %     matrix of gene names. That way, I can ask "where is pal-1 expressed
    %     in each strain". So, the final result will be a matrix where each row
    %     represents the position of a gene in each strain.

    
    %4. Create a matrix of gene names by their temporal order in each strain
    
    gene_order_each_strain = {};
    for i = 1:size(strains_in_3D,3)
        gene_order_each_strain(:,i) = Genes_no_lincs(selected_genes(gene_reorder(:,i)));
    end
    
    %5. Create a reference matrix of gene names (from the mean) arranged temporally
        
    gene_order_reference = Genes_no_lincs(selected_genes(gene_order_of_mean_expression_matrix));%Gene order in mean of 20 strains
    gene_names_order_distance_matrix = gene_order_reference;
    
%     try
%         load gene_pos_ind_in_means_matrix
%     catch
        gene_pos_ind_each_strain = [];
        for i = 1:size(gene_order_each_strain,2)
            for j = 1:size(gene_order_each_strain,1)
                gene_pos_ind_each_strain(j,i) = find(ismember(gene_order_each_strain(:,i),gene_order_reference(j)));
            end
        end
        
        save gene_pos_ind_in_means_matrix gene_pos_ind_each_strain     
%     end

          
    gene_pos_ind_ref = (1:size(gene_order_reference))';
    
    Distance1 = [];Distance2 = [];Distance = [];
    
    for i = 1:size(gene_pos_ind_each_strain,2)
        Distance1(:,i) = abs(gene_pos_ind_each_strain(:,i)- gene_pos_ind_ref);
        Distance2(:,i) = size(gene_pos_ind_ref,1)-abs(gene_pos_ind_each_strain(:,i)-gene_pos_ind_ref);
    end
    
    %%%The matrices 'Distance' and 'mean_distance' are already temporally
    %%%arranged. the arrangement is by the reference list used
    
    for j = 1:size(Distance1,2)
        for i = 1:size(Distance1,1)
            Distance(i,j) = min([Distance1(i,j) Distance2(i,j)]);
        end
    end
    
 
    mean_distance_temporal = mean(Distance,2);

    
    RW_Genes = {};
    Genes_sorted_strains = Genes_no_lincs(selected_genes(gene_order_of_mean_expression_matrix));
    
    ii = 1:step_size:size(mean_distance_temporal,1)-window_size;
    
    for k = 1:size(Distance,2)
        for i = 1:length(ii)
            for j = ii(1,i)
                RW_Genes(:,i) = Genes_sorted_strains(j:j+window_size);
            end
        end
    end
       
    save distance_RW_data mean_distance_temporal  Distance;
         
    save Gene_names_RW RW_Genes gene_names_order_distance_matrix gene_order_each_strain...
        gene_order_reference Genes_sorted_strains;

end

if (running_window_expression_pattern)
    
    
    selected_genes = find(isnan(meanCorr)==0);%genes with no NaNs in the 20 strains corr matrix
    
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
    
    Z_mean_expression_matrix = zscore(mean_expression_matrix(selected_genes,:),0,2);
    
    [pcg, zscoresg] = pca(Z_mean_expression_matrix);
    
    zavit = mod(atan2d(zscoresg(:,1),zscoresg(:,2)), 0);
    [~, gene_order_of_mean_expression_matrix] = sort(zavit);      
    
    correlation_with_start = corrcoef([ 1 0 0 0 0 0 0; Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:) ]');
    [m, i_m] = max(correlation_with_start(1,2:end));
    gene_order_of_mean_expression_matrix = [gene_order_of_mean_expression_matrix(i_m:end)', gene_order_of_mean_expression_matrix(1:i_m-1)'];%%%%
    [m, m_i]=max(Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:)');
    
    c = diag(corrcoef(m_i, 1:length(m_i)),1);
    
    if (c<0)
        gene_order_of_mean_expression_matrix = flipud(gene_order_of_mean_expression_matrix');
    end
    
    if (manual_move_genes_in_ZAVIT)        
       gene_order_of_mean_expression_matrix = circshift(gene_order_of_mean_expression_matrix',permutation_ind); 
    end
       
    mean_expression_matrix_raw = [];
    strt_num = 1:20:140;end_num = 20:20:140;
    for i = 1:7
        mean_expression_matrix_raw(:,i) = mean(values(:,strt_num(i):end_num(i)),2);
    end
    
    mean_expression_matrix_raw_ordered = mean_expression_matrix_raw(selected_genes(gene_order_of_mean_expression_matrix),:);
    
    Expression_ordered = Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:);
    
    Expression_ordered_log_values = mean_expression_matrix(selected_genes(gene_order_of_mean_expression_matrix),:);
        
    ii = 1:step_size:size(Z_mean_expression_matrix,1)-window_size;
    
    Running_window_expression = {};
    Running_window_expression_log_values = {};
    sum_Running_window_expression_log_values = [];
    mean_sum_Running_window_expression_log_values = [];
    
    for i = 1:length(ii)
        for j = ii(1,i)
           
            Running_window_expression{i} = Expression_ordered(j:j+window_size,:);
            Running_window_expression_log_values{i} = Expression_ordered_log_values(j:j+window_size,:);
            sum_Running_window_expression_log_values(:,i) = sum(Running_window_expression_log_values{i},2);
            mean_sum_Running_window_expression_log_values(i) = mean(sum(Running_window_expression_log_values{i},2));
        
        end
    end
    
    
    figure;
    for i = 1:size(Running_window_expression,2)
        subplot(5,11,i);
        imagesc(mean(Running_window_expression{i},1));colormap(cmap_color_blind2);
        
        set(gca,'xtick',[],'ytick',[])
    end
    
    figure;
    h = boxplot(sum_Running_window_expression_log_values,'colors','b','notch','on');
    set(h(7,:),'Visible','off')
    set(gca, 'FontSize', 10,'xticklabel',[])
    ylabel('Expression level');xlabel('Running window cluster');
    
    h2 = findobj(gca,'Tag','Box');
    for j=1:length(h2)
        patch(get(h2(j),'XData'),get(h2(j),'YData'),'r','FaceAlpha',.5);
    end;
            
    save running_window_expression_pattern Expression_ordered_log_values sum_Running_window_expression_log_values Running_window_expression Running_window_expression_log_values
    
end

if (replicates_analysis)
    
    %Here, we calculate the same as in the previous section, only for the 
    %four biological replicates
   
    selected_genes = find(isnan(meanCorr)==0); 
    
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
   
    try
        load MA_bioreplicate_norm
    catch
        
        biorepeat_UMI = real(-(4^5)*log(1-biorepeat_data_2/(4^5)));
        
        biorepeat_expression_tpm = [];
        for i = 1:size(biorepeat_data_2,2)
            biorepeat_expression_tpm(:,i) = biorepeat_UMI(:,i)*50000/(sum(biorepeat_UMI(:,i),1));
        end
        
        
        if (log_norm)
            biorepeat_expression_tpm_log = log2(1 + biorepeat_expression_tpm);
        else
            biorepeat_expression_tpm_log = biorepeat_expression_tpm;
        end
        
        biorepeat_expression_tpm_log = biorepeat_expression_tpm_log(1:20517,:);
        
        
        
        biorepeat_replicates_in_3D = [];
        for i = 1:4
            biorepeat_replicates_in_3D(:,:,i) =  biorepeat_expression_tpm_log(:,i:4:size(biorepeat_expression_tpm_log,2));
        end
        
        mean_expression_data_biorepeat = [];
        for i = 1:7
            mean_expression_data_biorepeat(:,i) = mean(biorepeat_replicates_in_3D(:,i,:),3);
        end
        
        save MA_bioreplicate_norm biorepeat_expression_tpm...
            biorepeat_expression_tpm_log biorepeat_replicates_in_3D...
            mean_expression_data_biorepeat
        
    end
    
    
        
    
    figure;imagesc(corr(biorepeat_expression_tpm_log(selected_genes,:)));
    colormap(cmap_color_blind2);colorbar
    set(gca,'xtick',[],'ytick',[]);

    
    Z_mean_expression_matrix_replicates = zscore(mean_expression_data_biorepeat(selected_genes,:),0,2);
    
    [pcg, zscoresg] = pca(Z_mean_expression_matrix_replicates);
    
    zavit = mod(atan2d(zscoresg(:,1),zscoresg(:,2)), 0);
    [~, gene_order_of_mean_expression_matrix_replicates] = sort(zavit);      
    
    correlation_with_start = corrcoef([ 1 0 0 0 0 0 0;...
        Z_mean_expression_matrix_replicates(gene_order_of_mean_expression_matrix_replicates,:) ]');
    [m, i_m] = max(correlation_with_start(1,2:end));
    gene_order_of_mean_expression_matrix_replicates = [gene_order_of_mean_expression_matrix_replicates(i_m:end)', gene_order_of_mean_expression_matrix_replicates(1:i_m-1)'];%%%%
    [m, m_i]=max(Z_mean_expression_matrix_replicates(gene_order_of_mean_expression_matrix_replicates,:)');
    
    c = diag(corrcoef(m_i, 1:length(m_i)),1);
    
    if (c<0)
        gene_order_of_mean_expression_matrix_replicates = flipud(gene_order_of_mean_expression_matrix_replicates');
    end
    
    if (manual_move_genes_in_ZAVIT)        
       gene_order_of_mean_expression_matrix_replicates = circshift(gene_order_of_mean_expression_matrix_replicates,permutation_ind); 
    end
     
    figure;imagesc(Z_mean_expression_matrix_replicates(gene_order_of_mean_expression_matrix_replicates,:));colormap(redblue);
    Genes_sorted = Genes_no_lincs(selected_genes(gene_order_of_mean_expression_matrix_replicates));
    
    
    ii = 1:step_size:size(Z_mean_expression_matrix_replicates,1)-window_size;
    RW_Genes_BR = {};
    for i = 1:length(ii)
        for j = ii(1,i)
            RW_Genes_BR(:,i) = Genes_sorted(j:j+window_size);
        end
    end
           
        each_replicate_fitted = [];PRS = [1 0 0 0 0 0 0];
        pcg = [];zscoresg = [];pcvarsg = [];zavit = [];gene_order = [];
        each_strain_fitted = [];zavit = [];
        
        for i = 1:size(biorepeat_replicates_in_3D,3);
            
            each_replicate_fitted(:,:,i) = zscore(biorepeat_replicates_in_3D(selected_genes,:,i),0,2);
            
            [pcg(:,:,i), zscoresg(:,:,i), pcvarsg(:,:,i)] = pca(each_replicate_fitted(:,:,i));
            zavit(:,i) = mod(atan2d(zscoresg(:,1,i),zscoresg(:,2,i)), 0) ;
            [~, gene_order(:,i)] = sort(zavit(:,i));
            
        end  
        
        %This produces a messy picture. To fix this, we initially find the gene
        %in each strain most correlated with the profile [1 0 0 0 0 0 0]
        ci = [];pi = [];exp_prof = [];
        
        for i = 1:size(each_replicate_fitted,3);
            [ci(:,i),pi(:,i)] = corr(each_replicate_fitted(:,:,i)',PRS(1,:)');
            [~,exp_prof(i)] = max(ci(:,i)');
        end
        
        each_replicate_fitted_and_sorted = [];gene_reorder = [];
        for i = 1:size(each_replicate_fitted,3)
         
            first_gene(i) = find(gene_order(:,i) == exp_prof(i));

            gene_reorder(:,i) = circshift(gene_order(:,i),(size(gene_order,1)-(first_gene(i)-1)));
            
            [~,m_i] = max(each_replicate_fitted(gene_reorder(:,i),:,i));
            c  =corrcoef(m_i, 1:length(m_i)); c = c(1,2);
            if (c<0.3)
                gene_reorder(:,i) = flipud(gene_reorder(:,i));
            end
            
            
        end
        
        if (manual_move_genes_in_ZAVIT)
            gene_reorder = circshift(gene_reorder,permutation_ind);
        end
        
        for i = 1:size(each_replicate_fitted,3)
            each_replicate_fitted_and_sorted(:,:,i) = each_replicate_fitted(gene_reorder(:,i),:,i);
        end
        
        figure;
        for i = 1:size(each_replicate_fitted_and_sorted,3)
            subplot(1,4,i);imagesc(each_replicate_fitted_and_sorted(:,:,i));
            colormap(redblue);
        end
        
        gene_order_each_replicate = {};
        for i = 1:size(gene_reorder,2);
            gene_order_each_replicate(:,i) = Genes_no_lincs(selected_genes(gene_reorder(:,i)));
        end
        
        gene_pos_ind = [];
        for i = 1:size(gene_order_each_replicate,2)
            for j = 1:size(gene_order_each_replicate,1)
                gene_pos_ind(j,i) = find(ismember(gene_order_each_replicate(:,i),Genes_sorted(j)));
            end
        end
        
        save gene_pos_ind_replicates gene_pos_ind;
                        
        gene_pos_ind_ref = (1:size(Genes_sorted))';
        
        Distance1 = [];Distance2 = [];Distance_BR = [];
        
        for i = 1:size(gene_pos_ind,2)
            Distance1(:,i) = abs(gene_pos_ind(:,i)- gene_pos_ind_ref);
            Distance2(:,i) = size(gene_pos_ind_ref,1)-abs(gene_pos_ind(:,i)-gene_pos_ind_ref);
        end
         
        
        %%%The matrices 'Distance' and 'mean_distance' are already temporally
        %%%arranged. the arrangement is by the reference list used
        
        for j = 1:size(gene_order_each_replicate,2)
            for i = 1:size(Distance1,1)
                Distance_BR(i,j) = min([Distance1(i,j) Distance2(i,j)]);
            end
        end
        
        mean_distance_temporal_br = mean(Distance_BR,2);          
                
    Expression_ordered_br = Z_mean_expression_matrix_replicates(...
        gene_order_of_mean_expression_matrix_replicates,:);
    
    Expression_ordered_log_values_br = mean_expression_data_biorepeat(...
        selected_genes(gene_order_of_mean_expression_matrix_replicates),:);
        
    ii = 1:step_size:size(Z_mean_expression_matrix_replicates,1)-window_size;
    
    Running_window_expression_br = {};
    Running_window_expression_log_values_br = {};
    sum_Running_window_expression_log_values_br = [];
    mean_sum_Running_window_expression_log_values_br = [];
    
    for i = 1:length(ii)
        for j = ii(1,i)
           
            Running_window_expression_br{i} = Expression_ordered_br(j:j+window_size,:);
            Running_window_expression_log_values_br{i} = Expression_ordered_log_values_br(j:j+window_size,:);
            sum_Running_window_expression_log_values_br(:,i) = sum(Running_window_expression_log_values_br{i},2);
            mean_sum_Running_window_expression_log_values_br(i) = mean(sum(Running_window_expression_log_values_br{i},2));
        end
    end


    figure;
    for i = 1:size(Running_window_expression_log_values_br,2)
        subplot(5,11,i);
        imagesc(mean(Running_window_expression_log_values_br{i},1));
        colormap(cmap_color_blind2);
    end
    
    figure;
    h = boxplot(sum_Running_window_expression_log_values_br,'colors','b','notch','on');
    set(h(7,:),'Visible','off')
    set(gca, 'FontSize', 10,'xticklabel',[])
    ylabel('Expression');xlabel('Running window cluster');
    
    h2 = findobj(gca,'Tag','Box');
    for j=1:length(h2)
        patch(get(h2(j),'XData'),get(h2(j),'YData'),'r','FaceAlpha',.5);
    end;
    
    save distance_RW_replicates_data mean_distance_temporal_RW_br mean_distance_temporal_br...
        gene_order_each_replicate Genes_sorted Expression_ordered_log_values_br;

end

if (distance_cleanup_and_running_window)
       
    randomize_order = 0;
   
    selected_genes = find(isnan(meanCorr)==0); 
    
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);

    mean_distance_temporal_noise = abs(mean_distance_temporal-mean_distance_temporal_br);
    mean_distance_temporal_clean = abs(mean_distance_temporal - mean_distance_temporal_noise);
    mean_distance_temporal_clean_br = abs(mean_distance_temporal_br - mean_distance_temporal_noise);

    if (randomize_order)
        rand_perm_ind  = randperm(length(selected_genes));
        mean_distance_temporal_clean = mean_distance_temporal_clean(rand_perm_ind);
    end
    
    ii = 1:step_size:size(mean_distance_temporal_clean,1)-window_size;
    mean_distance_temporal_RW_clean = [];
    
    for i = 1:length(ii)
        for j = ii(1,i)
            mean_distance_temporal_RW_clean(:,i) = mean_distance_temporal_clean(j:j+window_size);
        end
    end
    
    figure;
    h = boxplot(fliplr(mean_distance_temporal_RW_clean),'colors',[0 0 0],'notch','on','widths',0.4);
    axis([0 56 50 400])
    set(h([1:4 7],:),'Visible','off');set(h(5,:),'LineWidth',2)
    set(gca,'ytick',[],'yticklabel',[],'xtick',[],'xticklabel',[]);

    h2 = findobj(gca,'Tag','Box');

    for j=1:length(h2)
        patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.9 0.9 0.9],'FaceAlpha',.5);
    end;
    
    hold on
    plot(1:size(mean_distance_temporal_RW_clean,2),median(fliplr(mean_distance_temporal_RW_clean))...
        ,'LineWidth',2,'color','black');

    plot(0:56,repmat(mean(mean_distance_temporal_clean),1,57),'LineWidth',2,'color','black')
    oj1 = [];
    for i = 1:size(mean_distance_temporal_RW_clean,2)
        [~,oj1(i)] = kstest2(mean_distance_temporal_RW_clean(:,i),mean_distance_temporal_clean);
    end
    
    
    pj1 = find(mean(mean_distance_temporal_RW_clean)>mean(mean_distance_temporal_clean) & oj1<5e-10);
    pj2 = find(mean(mean_distance_temporal_RW_clean)<mean(mean_distance_temporal_clean) & oj1<5e-10);
     
    
    for i = 1:length(pj1);
        for j = pj1(i)
            patch(get(h2(j),'XData'),get(h2(j),'YData'),[1 0 0],'FaceAlpha',.5);
        end;
    end 
    
    for i = 1:length(pj2);
        for j = pj2(i)
            patch(get(h2(j),'XData'),get(h2(j),'YData'),[0 0 1],'FaceAlpha',.5);
        end;
    end  
        
    save clean_distance mean_distance_temporal_clean mean_distance_temporal_RW_clean

end

if (GO_analysis)
             
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
    
    
    gene_list_matrix_to_be_used = [TF_families signaling  GO_Matrix_new];
    gene_list_titles_matrix_to_be_used = [TF_family_names GS_names GO_Names_new'];
    
    %     remove redundant columns:
    
    z = triu(corr(gene_list_matrix_to_be_used),1);
    [a,b] = find(z>0.6);
    
    gene_list_matrix_to_be_used(:,b) = [];
    gene_list_titles_matrix_to_be_used(b) = [];
        
    GO_groups_expression = [];mean_distance_gene_lists = {};mean_mean_distance_gene_lists = [];    
    
    for i = 1:size(gene_list_matrix_to_be_used,2)
        
        GO_groups_expression(i,:) = sum(ismember(...
            RW_Genes,Genes_no_lincs(logical(gene_list_matrix_to_be_used(:,i)))));
        
        mean_distance_gene_lists{i} = mean_distance_temporal_clean(find(ismember...
            (gene_names_order_distance_matrix,Genes_no_lincs(logical(gene_list_matrix_to_be_used(:,i))))));
        
        mean_mean_distance_gene_lists(i) = mean(mean_distance_gene_lists{i});
        
    end
    
    % Remove GOs with less than 30 genes expressed in the dataset, and
    % filter for dynamic GOs
    x = sum(GO_groups_expression,2);
    x_thresh = 30;
    
    GO_groups_expression = GO_groups_expression(x>x_thresh,:);
    mask = genevarfilter(GO_groups_expression,'percentile',50);
    GO_groups_expression = GO_groups_expression(mask,:);
  
    mean_mean_distance_gene_lists = mean_mean_distance_gene_lists(x>x_thresh);
    mean_mean_distance_gene_lists = mean_mean_distance_gene_lists(mask);
    
    mean_distance_gene_lists = mean_distance_gene_lists(x>x_thresh);
    mean_distance_gene_lists = mean_distance_gene_lists(mask);
        
    gene_list_matrix_to_be_used = gene_list_matrix_to_be_used(:,x>x_thresh);
    gene_list_matrix_to_be_used = gene_list_matrix_to_be_used(:,mask);
    gene_list_titles_matrix_to_be_used = gene_list_titles_matrix_to_be_used(x>x_thresh);
    gene_list_titles_matrix_to_be_used = gene_list_titles_matrix_to_be_used(mask);
    
    Z = {};Z2 = [];
    for i = 1:size(gene_list_matrix_to_be_used,2)
        Z{i} = mean_expression_matrix(selected_genes(logical(gene_list_matrix_to_be_used(selected_genes,i))),:);
        Z2(i,:) = mean(Z{i},1);
    end
    
    Z2_s = zscore(Z2,0,2);
    GO_Names_new_ordered = gene_list_titles_matrix_to_be_used;
    GO_Names_new_ordered(isnan(mean(Z2_s,2))) = [];
    
    Z2_s(isnan(mean(Z2_s,2)),:) = [];
    [score,pc] = pca(Z2_s');
    zavit = mod(atan2d(score(:,1),score(:,2)), 0);
    [~, GO_order] = sort(zavit);
    
    GO_order = flipud(GO_order);
    
    GO_Names_new_ordered = GO_Names_new_ordered(GO_order);
    
    figure;imagesc(Z2_s(GO_order,:));colormap(cmap_color_blind2);colorbar;
    set(gca,'xtick',[],'ytick',1:length(GO_order),'yticklabel',GO_Names_new_ordered)
    
    figure;imagesc(corr(Z2_s(GO_order,:)'));colormap(cmap_color_blind2);colorbar;
    set(gca,'xtick',[],'ytick',[]);
       
    pv_GO = [];
    
    for i = 1:size(mean_distance_gene_lists,2)
        [~,pv_GO(i)] = kstest2(mean_distance_gene_lists{i},mean_distance_temporal_clean);
    end
    pv_GO = pv_GO(GO_order);
    
    ff = mean_mean_distance_gene_lists(GO_order);
    ff_unconserved = find(ff>mean(mean_distance_temporal_clean) & pv_GO < 0.05);
    ff_conserved = find(ff<mean(mean_distance_temporal_clean) & pv_GO < 0.05);
    
    figure;bar(ff,'y');axis tight;
    ylim([100 max(mean_mean_distance_gene_lists)]);
    hold on
    bar(ff_conserved,ff(ff_conserved),'b');
    bar(ff_unconserved,ff(ff_unconserved),'r');
%     set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    plot(0:length(ff)+1,repmat(mean(mean_distance_temporal_clean),1,length(ff)+2),'LineWidth',2,'color','black')
    
    ff_sig = [ff_unconserved ff_conserved];
    [fff,ffff] = sort(ff(ff_sig),'descend');
    ff_sig_sorted = ff_sig(ffff);
    pv_GO_sorted = pv_GO(ff_sig_sorted);
    figure;imagesc(-log10(pv_GO_sorted))
    
    
    figure;imagesc(Z2_s(GO_order(ff_sig_sorted),:));colormap(cmap_color_blind2);colorbar;
    set(gca,'xtick',[],'ytick',1:length(GO_order),'yticklabel',GO_Names_new_ordered(ff_sig_sorted))
    

    figure;
    bar(1:find(fff == min(ff(ff_unconserved))),fff(1:find(fff == min(ff(ff_unconserved)))),'r');
    hold on
    bar((find(fff == min(ff(ff_unconserved)))+1):length(fff),...
        fff((find(fff == min(ff(ff_unconserved)))+1):length(fff)),'b');
    ylim([min(mean_mean_distance_gene_lists)-10 max(mean_mean_distance_gene_lists)+10]);
    xlim([0 length(fff)+1]);
%     set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    hold on
    plot(0:length(fff)+1,repmat(mean(mean_distance_temporal_clean),1,length(fff)+2),'LineWidth',2,'color','black')

    

    
end

if (functional_groups_analysis) 
    
  
    selected_genes = find(isnan(meanCorr)==0);%genes with no NaNs in the 20 strains corr matrix
    
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
            
    gene_list_matrix_to_be_used = [TF_families(:,7) specific_lists(:,16) specific_lists(:,18)...
        ismember(Genes_no_lincs,missing_genes) specific_lists(:,17)];
    gene_list_titles_matrix_to_be_used = {'Homeodomain' 'Endoderm' 'Ectoderm' '"Missing Genes"' 'Mesoderm'};

    
    GO_groups_expression = [];
    
    for i = 1:size(gene_list_matrix_to_be_used,2)
        GO_groups_expression(i,:) = sum(ismember(RW_Genes,Genes_no_lincs(logical(gene_list_matrix_to_be_used(:,i)))));
    end
    x = sum(GO_groups_expression,2);
    
    p_thresh = 1;
    
    
    mean_distance_gene_lists = {};mean_mean_distance_gene_lists = [];
    std_mean_distance_gene_lists = [];
    size_mean_distance_gene_lists = [];grp = {};
    
    
    for i = 1:size(gene_list_matrix_to_be_used,2)
        mean_distance_gene_lists{i} = mean_distance_temporal_clean(find(ismember...
            (gene_names_order_distance_matrix,Genes_no_lincs(logical(gene_list_matrix_to_be_used(:,i))))));
        mean_mean_distance_gene_lists(i) = mean(mean_distance_gene_lists{i});
        std_mean_distance_gene_lists(i) = std(mean_distance_gene_lists{i});
        size_mean_distance_gene_lists(i) = size(mean_distance_gene_lists{i},1);
        grp{i} = repmat(i,size_mean_distance_gene_lists(i),1);
    end
    
    grpng = vertcat(grp{1:end});
    
    d = [];p = [];jojo = [];
    for i = 1:size(gene_list_matrix_to_be_used,2)
        jojo(i) = isempty(mean_distance_gene_lists{i});
        mean_distance_gene_lists(logical(jojo)) = {mean_distance_temporal_clean};
        [d(i),p(i)] = kstest2(mean_distance_temporal_clean,mean_distance_gene_lists{i});
    end
    
    mean_distance_gene_lists_all = vertcat(mean_distance_gene_lists{1:end});
    
    gene_list_significantly_unconserved = gene_list_titles_matrix_to_be_used(mean_mean_distance_gene_lists>mean(mean_distance_temporal_clean) & p<p_thresh);
    gene_list_significantly_conserved = gene_list_titles_matrix_to_be_used(mean_mean_distance_gene_lists<mean(mean_distance_temporal_clean) & p<p_thresh);
        
    %     draw only gene groups significantly more or less conserved:
    
    gene_list_significantly_unconserved_ind = find(ismember(gene_list_titles_matrix_to_be_used,gene_list_significantly_unconserved));
    gene_list_significantly_conserved_ind = find(ismember(gene_list_titles_matrix_to_be_used,gene_list_significantly_conserved));
    
    if (size(gene_list_significantly_unconserved_ind,1)>size(gene_list_significantly_unconserved_ind,2))
        gene_list_significantly_altered = [gene_list_significantly_unconserved_ind;gene_list_significantly_conserved_ind];
    else
        gene_list_significantly_altered = [gene_list_significantly_unconserved_ind gene_list_significantly_conserved_ind];
    end
    
    [cho,cho2] = sort(mean_mean_distance_gene_lists(gene_list_significantly_unconserved_ind),'descend');
    [ne_chievo,ne_chievo_2] = sort(mean_mean_distance_gene_lists(gene_list_significantly_conserved_ind),'descend');
    ne_huya = std_mean_distance_gene_lists(gene_list_significantly_unconserved_ind(cho2));
    pizdish = std_mean_distance_gene_lists(gene_list_significantly_conserved_ind(ne_chievo_2));
    
    tralala = gene_list_titles_matrix_to_be_used(gene_list_significantly_altered);
    
    
    figure;
    if (length(gene_list_significantly_unconserved_ind) >  0)
        bar(1:length(gene_list_significantly_unconserved_ind),mean_mean_distance_gene_lists(gene_list_significantly_unconserved_ind));
        hold on
        hArray = bar(length(gene_list_significantly_unconserved_ind)+1:length(gene_list_significantly_unconserved_ind)...
            +length(gene_list_significantly_conserved_ind),mean_mean_distance_gene_lists(gene_list_significantly_conserved_ind));
        set(hArray(1),'facecolor','red','edgecolor', 'red')
        %%%%%
        hold on
        errorbar([mean_mean_distance_gene_lists(gene_list_significantly_unconserved_ind)...
            mean_mean_distance_gene_lists(gene_list_significantly_conserved_ind)],...
            [std_mean_distance_gene_lists(gene_list_significantly_unconserved_ind)...
            std_mean_distance_gene_lists(gene_list_significantly_conserved_ind)],'xblack');
        
        
        %         axis tight
        set(gca,'xtick', 1:length(gene_list_significantly_altered), 'xticklabel', tralala, 'FontSize', 8)
        rotateXLabels(gca,45)
    else
        hArray = bar(length(gene_list_significantly_unconserved_ind)+1:length(gene_list_significantly_unconserved_ind)...
            +length(gene_list_significantly_conserved_ind),mean_mean_distance_gene_lists(gene_list_significantly_conserved_ind));
        set(hArray(1),'facecolor','blue','edgecolor', 'blue')
%         axis tight
        ylim([0 200])
        set(gca,'xtick', 1:length(gene_list_significantly_altered), 'xticklabel', tralala, 'FontSize', 8)
%         rotateXLabels(gca,45)
    end
    hold on;plot(0:length(tralala)+1,repmat(mean(mean_distance_temporal_clean),1,length(tralala)+2),'LineWidth',2,'color','black');

    figure;
    bar(1:5,mean_mean_distance_gene_lists([1 2 4 5 4]),'w');
    hold on    
    bar(6,mean(mean_distance_temporal_clean),'w')
    
    errorbar([mean_mean_distance_gene_lists([1 2 4 5 4]) mean(mean_distance_temporal_clean)],...
            [std_mean_distance_gene_lists([1 2 4 5 4]) std(mean_distance_temporal_clean)],'xblack');
        set(gca,'xtick',[],'ytick',[])
    
    [~,ne_znayu] = sort([mean_mean_distance_gene_lists(gene_list_significantly_unconserved_ind)...
        mean_mean_distance_gene_lists(gene_list_significantly_conserved_ind)],'descend');
    
    
    
    tralala = tralala(ne_znayu);
    ne_huya_pizdish = [ne_huya pizdish];
    ne_huya_pizdish = ne_huya_pizdish(ne_znayu);
    
    figure;
    if (length(gene_list_significantly_unconserved_ind) > 0)
        bar(1:length(gene_list_significantly_unconserved_ind),cho);
        hold on
        hArray = bar(length(gene_list_significantly_unconserved_ind)+1:length(gene_list_significantly_unconserved_ind)...
            +length(gene_list_significantly_conserved_ind),ne_chievo);
        set(hArray(1),'facecolor','red','edgecolor', 'red')
        
        hold on
        errorbar([cho ne_chievo],[ne_huya pizdish],'xblack');
        
        axis tight
        
        plot(0:length(tralala)+1,repmat(mean(mean_distance_temporal_clean),1,length(tralala)+2),'LineWidth',2,'color','black');
        set(gca,'xtick', 1:length(gene_list_significantly_altered), 'xticklabel', tralala, 'FontSize', 8)
        rotateXLabels(gca,45)
    else
        hArray = bar(length(gene_list_significantly_unconserved_ind)+1:length(gene_list_significantly_unconserved_ind)...
            +length(gene_list_significantly_conserved_ind),ne_chievo);
        set(hArray(1),'facecolor','red','edgecolor', 'red')
        axis tight
        set(gca,'xtick', 1:length(gene_list_significantly_altered), 'xticklabel', tralala, 'FontSize', 8)
        rotateXLabels(gca,45)
    end
      
    
end

if (figure_4_analysis)
    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
    
    missing_genes_matrix = ismember(Genes_no_lincs,missing_genes);
    
    Blast_tc_tpm_log = log2(1+Blast_tc_tpm);
    HZ_TH_corr = corr([mean_expression_matrix(selected_genes,:) Blast_tc_tpm_log(selected_genes,:)]);
    figure;
    imagesc(HZ_TH_corr(1:7,8:end));colorbar();
    set(gca,'xtick',[],'ytick',[])
    
    
    Germ_layer_exp_TH = [];
    
    Germ_layer_exp_TH(1,:) = mean(Blast_tc_tpm_log(selected_genes(specific_lists(selected_genes,16)),:));
    Germ_layer_exp_TH(2,:) = mean(Blast_tc_tpm_log(selected_genes(specific_lists(selected_genes,18)),:));    
    Germ_layer_exp_TH(3,:) = mean(Blast_tc_tpm_log(selected_genes(missing_genes_matrix(selected_genes,1)),:));
    Germ_layer_exp_TH(4,:) = mean(Blast_tc_tpm_log(selected_genes(specific_lists(selected_genes,17)),:));
    Germ_layer_exp_TH(5,:) = mean(Blast_tc_tpm_log(selected_genes(logical(TF_families(selected_genes,7)),:))); %(Homeodomain)
    
    Germ_layer_names = {'Endoderm' 'Ectoderm' 'Inegration Genes' 'Mesoderm' 'Homeodomain'};
    
    Germ_layer_exp_TH_stnd = [];
    
    for i = 1:size(Germ_layer_exp_TH,1);
        Germ_layer_exp_TH_stnd(i,:) = Germ_layer_exp_TH(i,:)./max(Germ_layer_exp_TH(i,:));
    end
    
    figure;
    imagesc(Germ_layer_exp_TH_stnd);
    set(gca,'xtick',[],'ytick',1:5,'yticklabel', Germ_layer_names)
    
    
    
    Germ_layer_exp_HZ = [];
    Germ_layer_exp_HZ(1,:) = mean(mean_expression_matrix(selected_genes(specific_lists(selected_genes,16)),:));
    Germ_layer_exp_HZ(2,:) = mean(mean_expression_matrix(selected_genes(specific_lists(selected_genes,18)),:));    
    Germ_layer_exp_HZ(3,:) = mean(mean_expression_matrix(selected_genes(missing_genes_matrix(selected_genes,1)),:));
    Germ_layer_exp_HZ(4,:) = mean(mean_expression_matrix(selected_genes(specific_lists(selected_genes,17)),:));
    Germ_layer_exp_HZ(5,:) = mean(mean_expression_matrix(selected_genes(logical(TF_families(selected_genes,7)),:)));

    Germ_layer_exp_HZ_stnd = [];
    
    for i = 1:size(Germ_layer_exp_HZ,1);
        Germ_layer_exp_HZ_stnd(i,:) = Germ_layer_exp_HZ(i,:)./max(Germ_layer_exp_HZ(i,:));
    end
    
    figure;imagesc(Germ_layer_exp_HZ_stnd);set(gca,'xtick',[],'ytick',1:5,'yticklabel', Germ_layer_names)
           
    Endoderm_ind_RW = {};Ectoderm_ind_RW = {};Mesoderm_ind_RW = {};Missing_ind_RW = {};
    for i = 1:size(RW_Genes,2)
        
        Endoderm_ind_RW{i} = find(ismember(RW_Genes(:,i), Genes_no_lincs(specific_lists(:,16))));
        Ectoderm_ind_RW{i} = find(ismember(RW_Genes(:,i), Genes_no_lincs(specific_lists(:,18))));
        Mesoderm_ind_RW{i} = find(ismember(RW_Genes(:,i), Genes_no_lincs(specific_lists(:,17))));
        Missing_ind_RW{i} = find(ismember(RW_Genes(:,i), Genes_no_lincs(missing_genes_matrix(:,1))));
        HD_ind_RW{i} = find(ismember(RW_Genes(:,i), Genes_no_lincs(logical(TF_families(:,7)))));
        
    end
    
    Endoderm_exp_RW = {};Ectoderm_exp_RW = {};Mesoderm_exp_RW = {};Missing_exp_RW = {};HD_exp_RW = {};
    for i = 1:size(RW_Genes,2)
        
        Endoderm_exp_RW{i} = sum_Running_window_expression_log_values(Endoderm_ind_RW{i},i);
        Ectoderm_exp_RW{i} = sum_Running_window_expression_log_values(Ectoderm_ind_RW{i},i);
        Mesoderm_exp_RW{i} = sum_Running_window_expression_log_values(Mesoderm_ind_RW{i},i);
        Missing_exp_RW{i} = sum_Running_window_expression_log_values(Missing_ind_RW{i},i);
        HD_exp_RW{i} = sum_Running_window_expression_log_values(HD_ind_RW{i},i);
        
        mean_Endoderm_exp_RW(i) = sum(Endoderm_exp_RW{i});
        mean_Ectoderm_exp_RW(i) = sum(Ectoderm_exp_RW{i});
        mean_Mesoderm_exp_RW(i) = sum(Mesoderm_exp_RW{i});
        mean_Missing_exp_RW(i) = sum(Missing_exp_RW{i});
        mean_HD_exp_RW(i) = sum(HD_exp_RW{i});
        
    end
    
    mean_Endoderm_exp_RW_stnd = mean_Endoderm_exp_RW./max(mean_Endoderm_exp_RW);
    mean_Ectoderm_exp_RW_stnd = mean_Ectoderm_exp_RW./max(mean_Ectoderm_exp_RW);
    mean_Mesoderm_exp_RW_stnd = mean_Mesoderm_exp_RW./max(mean_Mesoderm_exp_RW);
    mean_Missing_exp_RW_stnd = mean_Missing_exp_RW./max(mean_Missing_exp_RW);
    mean_HD_exp_RW_stnd = mean_HD_exp_RW./max(mean_HD_exp_RW);
    
    mean_germ_and_such_exp_stnd = [mean_Endoderm_exp_RW_stnd;...
        mean_Ectoderm_exp_RW_stnd;mean_Mesoderm_exp_RW_stnd;...
        mean_Missing_exp_RW_stnd;mean_HD_exp_RW_stnd;];
        
    
    figure;
    for i = 1:size(mean_germ_and_such_exp_stnd,1)
        subplot(5,1,i)
        imagesc(mean_germ_and_such_exp_stnd(i,:));
        colormap(cmap_color_blind2);colorbar;
        set(gca,'xtick',[],'ytick',1:5,'yticklabel', Germ_layer_names(i))
    end

 
end

if (control_analysis)
    
    longer_analysis             = 1;
    longer_analysis_replicates  = 1;
    longer_analysis_cleanup     = 1;
    
    manual_move_genes_in_ZAVIT  = 1;
    permutation_ind = -925;
    step_size = 20;window_size = 40;
    
    PRS = [1 0 0 0 0 0 0];
    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
    
    Z_mean_expression_matrix = zscore(mean_expression_matrix(selected_genes,:),0,2);
    
    try
        load kclusters
    catch
        [kclust_ind,C] = kmeans(Z_mean_expression_matrix,10);
        
        save kclusters kclust_ind C
    end
    
    
    [pcK,scoreK] = pca(C);
    zavitK = mod(atan2d(scoreK(:,1),scoreK(:,2)), 0) ;
    [~, K_order] = sort(zavitK);K_order = flipud(K_order);
    
    figure;
    for i = 1:size(C,1)
        subplot(2,5,i)
        plot(C(K_order(i),:),'LineWidth',2,'color','black');
        axis tight
        set(gca,'xtick',1:7,'xticklabel',[],'yticklabel',[])
    end
    
    kmeans_clusters_gene_names = {};kmeans_clusters_gene_ind_in_dist_matrix = {};
    kmeans_clusters_distance = {};kmeans_clusters_expression = {};permutation_ind_kmeans = {};
    kmeans_extracted_distance = {};kmeans_extracted_expression = {};
    
    
    
    for i = 1:max(kclust_ind)
        kmeans_clusters_gene_names{i} = Genes_no_lincs(selected_genes(kclust_ind == i));
        
        kmeans_clusters_gene_ind_in_dist_matrix{i} = ...
            find(ismember(gene_names_order_distance_matrix,kmeans_clusters_gene_names{i}));
        kmeans_clusters_expression{i} = Expression_ordered_log_values(kmeans_clusters_gene_ind_in_dist_matrix{i},:);
        
        
        try
            load permutation_ind_kmeans_2
        catch
            permutation_ind_kmeans{i} = randperm(size(kmeans_clusters_expression{i},1));
            
            save permutation_ind_kmeans_2 permutation_ind_kmeans
        end
        
        kmeans_clusters_gene_names_2{i} = kmeans_clusters_gene_names{i}(permutation_ind_kmeans{i}(1:100),:);
        kmeans_extracted_expression{i} = kmeans_clusters_expression{i}(permutation_ind_kmeans{i}(1:100),:);
        
    end
    
    kmeans_clusters_gene_names_all = vertcat(kmeans_clusters_gene_names_2{1:end});
    kmeans_extracted_expression_all = vertcat(kmeans_extracted_expression{1:end});
    
    
    if (longer_analysis)
        
        selected_genes = [];
        selected_genes = find(ismember(Genes_no_lincs,kmeans_clusters_gene_names_all));
        
        %1. ZAVIT for each strain seperately
        pcg = []; zscoresg = [];zavit = [];gene_order = [];
        for i = 1:size(strains_in_3D,3);
            
            each_strain_fitted(:,:,i) = zscore(strains_in_3D(selected_genes,:,i),0,2);
            
            [pcg(:,:,i), zscoresg(:,:,i)] = pca(each_strain_fitted(:,:,i));
            zavit(:,i) = mod(atan2d(zscoresg(:,1,i),zscoresg(:,2,i)), 0) ;
            [~, gene_order(:,i)] = sort(zavit(:,i));
            
        end
        
        %This produces a messy picture. To fix this, we initially find the gene
        %in each strain most correlated with the profile [1 0 0 0 0 0 0]
        
        ci = [];pi = [];exp_prof = [];
        
        for i = 1:size(each_strain_fitted,3);
            [ci(:,i),pi(:,i)] = corr(each_strain_fitted(:,:,i)',PRS(1,:)');
            [~,exp_prof(i)] = max(ci(:,i)');
        end
        %
        each_strain_fitted_and_sorted = [];gene_reorder = [];
        for i = 1:size(each_strain_fitted,3)
            
            first_gene(i) = find(gene_order(:,i) == exp_prof(i));
            
            gene_reorder(:,i) = circshift(gene_order(:,i),(size(gene_order,1)-(first_gene(i)-1)));
            
            [~,m_i] = max(each_strain_fitted(gene_reorder(:,i),:,i));
            c  =corrcoef(m_i, 1:length(m_i)); c = c(1,2);
            if (c<0.3)
                gene_reorder(:,i) = flipud(gene_reorder(:,i));
            end
            
            
        end
        
        if (manual_move_genes_in_ZAVIT)
            gene_reorder = circshift(gene_reorder,permutation_ind);
        end
        
        
        for i = 1:size(strains_in_3D,3);
            each_strain_fitted_and_sorted(:,:,i) = each_strain_fitted(gene_reorder(:,i),:,i);
        end
        
        
        figure;
        for i  = 1:size(each_strain_fitted_and_sorted,3)
            subplot(4,5,i)
            imagesc(each_strain_fitted_and_sorted(:,:,i));colormap(redblue)
        end
        
        %2. ZAVIT for the means matrix
        
        Z_mean_expression_matrix = zscore(mean_expression_matrix(selected_genes,:),0,2);
        [pcg, zscoresg] = pca(Z_mean_expression_matrix);
        
        zavit = mod(atan2d(zscoresg(:,1),zscoresg(:,2)), 0);
        [~, gene_order_of_mean_expression_matrix] = sort(zavit);
        
        
        correlation_with_start = corrcoef([ 1 0 0 0 0 0 0; Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:) ]');
        [m, i_m] = max(correlation_with_start(1,2:end));
        gene_order_of_mean_expression_matrix = [gene_order_of_mean_expression_matrix(i_m:end)', gene_order_of_mean_expression_matrix(1:i_m-1)'];%%%%
        [m, m_i]=max(Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:)');
        
        c = diag(corrcoef(m_i, 1:length(m_i)),1);
        
        if (c<0)
            gene_order_of_mean_expression_matrix = flipud(gene_order_of_mean_expression_matrix');
        end
        
        if (manual_move_genes_in_ZAVIT)
            gene_order_of_mean_expression_matrix = circshift(gene_order_of_mean_expression_matrix,permutation_ind);
        end
        
        figure;imagesc(Z_mean_expression_matrix(gene_order_of_mean_expression_matrix,:));colormap(redblue);colorbar;
        set(gca,'ytick',[200 400 600 800 1000],'yticklabel',[],'xtick',[])
        
        gene_order_each_strain = {};
        for i = 1:size(strains_in_3D,3)
            gene_order_each_strain(:,i) = Genes_no_lincs(selected_genes(gene_reorder(:,i)));
        end
        
        %3B. Create a reference matrix of gene names (from the mean) arranged temporally
        
        gene_order_reference = Genes_no_lincs(selected_genes(gene_order_of_mean_expression_matrix));%Gene order in mean of 20 strains
        
        gene_pos_ind_each_strain = [];
        for i = 1:size(gene_order_each_strain,2)
            for j = 1:size(gene_order_each_strain,1)
                gene_pos_ind_each_strain(j,i) = find(ismember(gene_order_each_strain(:,i),gene_order_reference(j)));
            end
        end
        gene_pos_ind_ref = (1:size(gene_order_reference))';
        
        Distance1 = [];Distance2 = [];Distance = [];
        
        for i = 1:size(gene_pos_ind_each_strain,2)
            Distance1(:,i) = abs(gene_pos_ind_each_strain(:,i)- gene_pos_ind_ref);
            Distance2(:,i) = size(gene_pos_ind_ref,1)-abs(gene_pos_ind_each_strain(:,i)-gene_pos_ind_ref);
        end
        
        %%%The matrices 'Distance' and 'mean_distance' are already temporally
        %%%arranged. the arrangement is by the reference list used
        
        for j = 1:size(Distance1,2)
            for i = 1:size(Distance1,1)
                Distance(i,j) = min([Distance1(i,j) Distance2(i,j)]);
            end
        end
        
        
        mean_distance_temporal = mean(Distance,2);
        
        
        ii = 1:step_size:size(kmeans_clusters_gene_names_all,1)-window_size;
        RW_Genes = {};
        Genes_sorted_strains = Genes_no_lincs(selected_genes(gene_order_of_mean_expression_matrix));
        
        for i = 1:length(ii)
            for j = ii(1,i)
                RW_Genes(:,i) = Genes_sorted_strains(j:j+window_size);
            end
        end
        
        
    end
    
    if (longer_analysis_replicates)
        
        selected_genes = [];
        selected_genes = find(ismember(Genes_no_lincs,kmeans_clusters_gene_names_all));
        
        %1. ZAVIT for each strain seperately
        
        
        gene_order_rep = [];zavit = [];pcg = [];zscoresg = [];
        for i = 1:size(biorepeat_replicates_in_3D,3);
            
            each_rep_fitted(:,:,i) = zscore(biorepeat_replicates_in_3D(selected_genes,:,i),0,2);
            
            [pcg(:,:,i), zscoresg(:,:,i)] = pca(each_rep_fitted(:,:,i));
            zavit(:,i) = mod(atan2d(zscoresg(:,1,i),zscoresg(:,2,i)), 0) ;
            [~, gene_order_rep(:,i)] = sort(zavit(:,i));
            
        end
        
        %This produces a messy picture. To fix this, we initially find the gene
        %in each strain most correlated with the profile [1 0 0 0 0 0 0]
        
        ci = [];pi = [];exp_prof = [];
        
        for i = 1:size(each_rep_fitted,3);
            [ci(:,i),pi(:,i)] = corr(each_rep_fitted(:,:,i)',PRS(1,:)');
            [~,exp_prof(i)] = max(ci(:,i)');
        end
        
        each_rep_fitted_and_sorted = [];gene_reorder_rep = [];
        
        for i = 1:size(each_rep_fitted,3)
            
            first_gene(i) = find(gene_order_rep(:,i) == exp_prof(i));
            
            gene_reorder_rep(:,i) = circshift(gene_order_rep(:,i),(size(gene_order_rep,1)-(first_gene(i)-1)));
            
            [~,m_i] = max(each_rep_fitted(gene_reorder_rep(:,i),:,i));
            c  =corrcoef(m_i, 1:length(m_i)); c = c(1,2);
            if (c<0.3)
                gene_reorder_rep(:,i) = flipud(gene_reorder_rep(:,i));
            end
        end
        
        if (manual_move_genes_in_ZAVIT)
            gene_reorder_rep = circshift(gene_reorder_rep,permutation_ind);
        end
        
        for i = 1:size(each_rep_fitted,3);
            each_rep_fitted_and_sorted(:,:,i) = each_rep_fitted(gene_reorder_rep(:,i),:,i);
        end
        
        
        figure;
        for i  = 1:4
            subplot(1,4,i)
            imagesc(each_rep_fitted_and_sorted(:,:,i));colormap(redblue)
        end
        
        %2. ZAVIT for the means matrix
        
        Z_mean_expression_matrix_br = zscore(mean_expression_data_biorepeat(selected_genes,:),0,2);
        [pcg, zscoresg] = pca(Z_mean_expression_matrix_br);
        
        zavit = mod(atan2d(zscoresg(:,1),zscoresg(:,2)), 0);
        [~, gene_order_of_mean_expression_matrix_br] = sort(zavit);
        
        
        correlation_with_start = corrcoef([ 1 0 0 0 0 0 0; Z_mean_expression_matrix(gene_order_of_mean_expression_matrix_br,:) ]');
        [m, i_m] = max(correlation_with_start(1,2:end));
        gene_order_of_mean_expression_matrix_br = [gene_order_of_mean_expression_matrix_br(i_m:end)', gene_order_of_mean_expression_matrix_br(1:i_m-1)'];%%%%
        [m, m_i]=max(Z_mean_expression_matrix_br(gene_order_of_mean_expression_matrix_br,:)');
        
        c = diag(corrcoef(m_i, 1:length(m_i)),1);
        
        if (c<0)
            gene_order_of_mean_expression_matrix_br = flipud(gene_order_of_mean_expression_matrix_br');
        end
        
        if (manual_move_genes_in_ZAVIT)
            gene_order_of_mean_expression_matrix_br = circshift(gene_order_of_mean_expression_matrix_br,permutation_ind);
        end
        
        figure;imagesc(Z_mean_expression_matrix_br(gene_order_of_mean_expression_matrix_br,:));colormap(redblue);colorbar;
        set(gca,'ytick',[200 400 600 800 1000],'yticklabel',[],'xtick',[])
        
        gene_order_each_strain_br = {};
        for i = 1:size(each_rep_fitted,3)
            gene_order_each_strain_br(:,i) = Genes_no_lincs(selected_genes(gene_reorder_rep(:,i)));
        end
        
        %3B. Create a reference matrix of gene names (from the mean) arranged temporally
        
        gene_order_reference_br = Genes_no_lincs(selected_genes(gene_order_of_mean_expression_matrix_br));%Gene order in mean of 20 strains
        
        gene_pos_ind_each_rep = [];
        for i = 1:size(gene_order_each_strain_br,2)
            for j = 1:size(gene_order_each_strain_br,1)
                gene_pos_ind_each_rep(j,i) = find(ismember(gene_order_each_strain_br(:,i),gene_order_reference_br(j)));
            end
        end
        gene_pos_ind_ref = (1:size(gene_order_reference_br))';
        
        Distance1 = [];Distance2 = [];Distance_br = [];
        
        for i = 1:size(gene_pos_ind_each_rep,2)
            Distance1(:,i) = abs(gene_pos_ind_each_rep(:,i)- gene_pos_ind_ref);
            Distance2(:,i) = size(gene_pos_ind_ref,1)-abs(gene_pos_ind_each_rep(:,i)-gene_pos_ind_ref);
        end
        
        %%%The matrices 'Distance' and 'mean_distance' are already temporally
        %%%arranged. the arrangement is by the reference list used
        
        for j = 1:size(Distance1,2)
            for i = 1:size(Distance1,1)
                Distance_br(i,j) = min([Distance1(i,j) Distance2(i,j)]);
            end
        end
        
        
        mean_distance_temporal_br = mean(Distance_br,2);
        
        
        
        
    end
    
    if (longer_analysis_cleanup)
        
        
        mean_distance_temporal_noise = abs(mean_distance_temporal - mean_distance_temporal_br);
        mean_distance_temporal_clean = abs(mean_distance_temporal - mean_distance_temporal_noise);
        mean_distance_temporal_clean_br = abs(mean_distance_temporal_br - mean_distance_temporal_noise);
        
        
        ii = 1:step_size:size(mean_distance_temporal_clean,1)-window_size;
        mean_distance_temporal_RW_clean = [];
        
        for i = 1:length(ii)
            for j = ii(1,i)
                mean_distance_temporal_RW_clean(:,i) = mean_distance_temporal_clean(j:j+window_size);
            end
        end
        
        figure;
        h = boxplot(fliplr(mean_distance_temporal_RW_clean),'colors',[0 0 0],'notch','on','widths',0.4);
        axis([0 size(mean_distance_temporal_RW_clean,2)+1 0 100])
        set(h([1:4 7],:),'Visible','off');set(h(5,:),'LineWidth',2)
        set(gca,'ytick',[],'yticklabel',[],'xtick',[],'xticklabel',[]);
        
        h2 = findobj(gca,'Tag','Box');
        
        for j=1:length(h2)
            patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.9 0.9 0.9],'FaceAlpha',.5);
        end;
        
        hold on
        plot(1:size(mean_distance_temporal_RW_clean,2),median(fliplr(mean_distance_temporal_RW_clean))...
            ,'LineWidth',2,'color','black');
        
        plot(0:size(mean_distance_temporal_RW_clean,2)+1,repmat(mean(mean_distance_temporal_clean)...
            ,1,size(mean_distance_temporal_RW_clean,2)+2),'LineWidth',2,'color','black')
        
        
        oj1 = [];
        for i = 1:size(mean_distance_temporal_RW_clean,2)
            [~,oj1(i)] = kstest2(mean_distance_temporal_RW_clean(:,i),mean_distance_temporal_clean);
        end
        
        pthr = 0.0005;
        
        pj1 = find(mean(mean_distance_temporal_RW_clean)>mean(mean_distance_temporal_clean) & oj1<pthr);
        pj2 = find(mean(mean_distance_temporal_RW_clean)<mean(mean_distance_temporal_clean) & oj1<pthr);
        
        
        for i = 1:length(pj1);
            for j = pj1(i)
                patch(get(h2(j),'XData'),get(h2(j),'YData'),[1 0 0],'FaceAlpha',.5);
            end;
        end
        
        for i = 1:length(pj2);
            for j = pj2(i)
                patch(get(h2(j),'XData'),get(h2(j),'YData'),[0 0 1],'FaceAlpha',.5);
            end;
        end
        
        
        %%%%%replicates:
        mean_distance_temporal_RW_clean_cont = [];
        for i = 1:length(ii)
            for j = ii(1,i)
                mean_distance_temporal_RW_clean_cont(:,i) = mean_distance_temporal_clean_br(j:j+window_size);
            end
        end
        
        figure;
        h = boxplot(fliplr(mean_distance_temporal_RW_clean_cont),'colors',[0 0 0],'notch','on','widths',0.4);
        axis([0 size(mean_distance_temporal_RW_clean_cont,2)+1 0 100])
        set(h([1:4 7],:),'Visible','off');set(h(5,:),'LineWidth',2)
        set(gca,'ytick',[],'yticklabel',[],'xtick',[],'xticklabel',[]);
        
        h2 = findobj(gca,'Tag','Box');
        
        for j=1:length(h2)
            patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.9 0.9 0.9],'FaceAlpha',.5);
        end;
        
        hold on
        plot(1:size(mean_distance_temporal_RW_clean_cont,2),median(fliplr(mean_distance_temporal_RW_clean_cont))...
            ,'LineWidth',2,'color','black');
        
        plot(0:size(mean_distance_temporal_RW_clean_cont,2)+1,repmat(mean(mean_distance_temporal_clean_br)...
            ,1,size(mean_distance_temporal_RW_clean_cont,2)+2),'LineWidth',2,'color','black')
        
        
        oj1 = [];
        for i = 1:size(mean_distance_temporal_RW_clean_cont,2)
            [~,oj1(i)] = kstest2(mean_distance_temporal_RW_clean_cont(:,i),mean_distance_temporal_clean_br);
        end
        
        
        pj1 = find(mean(mean_distance_temporal_RW_clean)>mean(mean_distance_temporal_clean_br) & oj1<pthr);
        pj2 = find(mean(mean_distance_temporal_RW_clean)<mean(mean_distance_temporal_clean_br) & oj1<pthr);
        
        
        for i = 1:length(pj1);
            for j = pj1(i)
                patch(get(h2(j),'XData'),get(h2(j),'YData'),[1 0 0],'FaceAlpha',.5);
            end;
        end
        
        for i = 1:length(pj2);
            for j = pj2(i)
                patch(get(h2(j),'XData'),get(h2(j),'YData'),[0 0 1],'FaceAlpha',.5);
            end;
        end
        
    end
    
    
end

if (GO_enrichment_for_gene_subgroups)
    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
        
    gene_list_matrix_to_be_used = [TF_families signaling  GO_Matrix_new];
    gene_list_titles_matrix_to_be_used = [TF_family_names GS_names GO_Names_new'];
    
    my_subgroup = missing_genes;
    
    GO_Gene_names = {};enrichment_sum_for_subgroup = [];GO_yn = [];GO_pv = [];
    for i = 1:size(gene_list_matrix_to_be_used,2);
        GO_Gene_names{i} = Genes_no_lincs(selected_genes(logical(...
            gene_list_matrix_to_be_used(selected_genes,i))));
        
        enrichment_sum_for_subgroup(i) = sum(ismember(GO_Gene_names{i},my_subgroup));
        
        GO_pv(i) = hygepdf(enrichment_sum_for_subgroup(i),...
            length(selected_genes),size(GO_Gene_names{i},1),length(my_subgroup));

    end
    
    enriched_GO_Names = gene_list_titles_matrix_to_be_used(GO_pv<0.05);
    
    figure;imagesc(-log10(GO_pv(GO_pv<0.05)));
    set(gca,'xtick',1:sum(GO_pv<0.05),'xticklabel',enriched_GO_Names);
    rotateXLabels(gca(),-90)
    
end

if (RW_to_stage_sorting)

    
    selected_genes = find(isnan(meanCorr)==0); %the genes with no NaNs in the correlation
    mask = genevarfilter(values_tpm_log(selected_genes,:), 'percentile', 30);
    selected_genes = selected_genes(mask == 1);
    
    a_sort = [];b_sort = [];Top_Genes = {};
    for i = 1:size(mean_expression_matrix,2)
        [a_sort(:,i),b_sort(:,i)] = sort(mean_expression_matrix(:,i),'descend');
        Top_Genes(:,i) = Genes_no_lincs(b_sort(1:1000,i));
    end
    
    sum_of_top_genes_in_RW = [];
    for j = 1:size(Top_Genes,2)
        for i = 1:size(RW_Genes,2)
            sum_of_top_genes_in_RW(j,i) = sum(ismember(RW_Genes(:,i), Top_Genes(:,j)));
        end
    end
    
    
    sum_of_top_genes_in_RW_stnd = [];
    for i = 1:size(sum_of_top_genes_in_RW,2)
        sum_of_top_genes_in_RW_stnd(:,i) = sum_of_top_genes_in_RW(:,i)/sum(sum_of_top_genes_in_RW(:,i));
    end
    
        
    [a_max,b_max] = max(sum_of_top_genes_in_RW_stnd,[],1);
   
    figure;imagesc(fliplr(b_max));colormap(paruly);
    set(gca,'xtick',[],'ytick',[]);
    
end


