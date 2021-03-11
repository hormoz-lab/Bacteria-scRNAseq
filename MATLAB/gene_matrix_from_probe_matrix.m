function [counts, dat, CBs, gene_name] = gene_matrix_from_probe_matrix(input_path, outdir, map_method)
%gene_matrix_from_probe_matrix generates a single-cell count matrix of
%genes from a single-cell count matrix of probes
%
%   gene_matrix_from_probe_matrix(INPUT_PATH, OUTDIR, MAP_METHOD) uses the
%   filtered_feature_bc_matrix.h5 and molecule_info.h5 files in INPUT_PATH 
%   for analysis and saves the result in .MTX format to OUTDIR. The number 
%   of columns in the input single-cell matrix file is equal to the number 
%   of probes. The number of columns in the output single-cell matrix file 
%   is equal to the number of genes. Set MAP_METHOD='max_probe_in_cell' to
%   use the maximally expressed probe for each gene and cell. Set
%   MAP_METHOD='median_probe_in_bulk' to use the same probe for each gene
%   across all cells. 
%
%   Author: Duluxan Sritharan. Hormoz Lab. Harvard Medical School.

    assert(~(exist(outdir, 'dir')==7), 'Output path already exists');
    
    if (~ismember(map_method, {'max_probe_in_cell'; 'median_probe_in_bulk'}))
        error('Unrecognized map_method: %s', map_method);
    end

    [counts, CBs, probe_name] = load_10x_h5_matrix(sprintf('%s/filtered_feature_bc_matrix.h5', input_path));    
    fprintf('Analysis for %d common cells using method=%s\n', length(CBs), map_method);

    dat   = parse_molecule_info(sprintf('%s/molecule_info.h5', input_path), CBs, length(probe_name)); 
    assert(isequal(dat.umi_counts, full(counts)));
    clear counts;
    
    split_loc = cellfun(@(x) strfind(x, '_'), probe_name, 'un', false);
    [gene_name, probe_id] = cellfun(@(x,s) deal(x(1:s(end)-1), x(s(end)+1:end)), probe_name, split_loc, 'un', false);
    probe_id = cellfun(@str2num, probe_id);
    [gene_name, ~, which_gene] = unique(gene_name, 'stable');
    
    reads_per_cell  = sum(dat.read_counts,2);
    umis_per_cell   = sum(dat.umi_counts,2);
    probes_per_cell = sum(logical(dat.umi_counts),2);
    genes_per_cell  = sum(logical(accumarray([dat.barcode_id, which_gene(dat.gene_id)], 1, [], [], [], true)),2);
    
    sum_probe = arrayfun(@(i) sum(dat.umi_counts(:,which_gene==i), 2), [1:length(gene_name)], 'un', false);
    sum_probe = horzcat(sum_probe{:});
    
    if (isequal(map_method, 'max_probe_in_cell'))
        [max_probe_umis, which_probe_max] = arrayfun(@(i) max(dat.umi_counts(:,which_gene==i), [], 2), [1:length(gene_name)], 'un', false);
        max_probe_reads = cellfun(@(g,w) g(sub2ind(size(g), [1:size(g,1)]', w)), arrayfun(@(i) dat.read_counts(:,which_gene==i), [1:length(gene_name)], 'un', false), which_probe_max, 'un', false);
        max_probe_umis = horzcat(max_probe_umis{:});
        which_probe_max = horzcat(which_probe_max{:});
        max_probe_reads = horzcat(max_probe_reads{:});
    else    
        bulk_gene = arrayfun(@(i) sum(dat.umi_counts(:,which_gene==i), 1), [1:length(gene_name)], 'un', false);
        which_probe_bulk_median = cellfun(@(x,y) y(find(abs(x-median(x))==min(abs(x-median(x))), 1, 'first')), bulk_gene, ...
                                          arrayfun(@(i) find(which_gene==i), [1:length(gene_name)], 'un', false));
        bulk_median_probe_reads = dat.read_counts(:, which_probe_bulk_median);
    end

    counts.by_probe      = dat.umi_counts;
    counts.by_gene_total = sum_probe;
    
    if (isequal(map_method, 'max_probe_in_cell'))    
        counts.by_gene_max_probe = max_probe_umis;
        counts.which_probe_max = which_probe_max;
        counts.max_probe_reads = max_probe_reads;
    else
        counts.by_median_bulk_probe = dat.umi_counts(:,which_probe_bulk_median);
        counts.which_probe_bulk_median = which_probe_bulk_median;
        counts.bulk_median_reads = bulk_median_probe_reads;
    end    
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,4,1);
    semilogy(sort(genes_per_cell, 'descend'));
    axis tight;
    ylim([1 max(genes_per_cell)])
    xlabel('Cell (ranked');
    ylabel('Unique Genes');
    
    subplot(2,4,2);
    semilogy(sort(probes_per_cell, 'descend'));
    axis tight;
    ylim([1 max(probes_per_cell)])
    xlabel('Cell (ranked');
    ylabel('Unique Probes');    
      
    subplot(2,4,3);
    if (isequal(map_method, 'max_probe_in_cell'))    
        semilogy([sort(umis_per_cell, 'descend') sort(sum(counts.by_gene_max_probe,2), 'descend')]);
        axis tight;
        ylim([1 max(umis_per_cell)])
        xlabel('Cell (ranked');
        ylabel('UMIs');
        legend({'All Probes'; 'Max Probe'});
        title(sprintf('%% UMIs Retained=%5.2f', full(sum(counts.by_gene_max_probe(:))/sum(umis_per_cell))*100));
    else
        semilogy([sort(umis_per_cell, 'descend') sort(sum(counts.by_median_bulk_probe,2), 'descend')]);
        axis tight;
        ylim([1 max(umis_per_cell)])
        xlabel('Cell (ranked');
        ylabel('UMIs');
        legend({'All Probes'; 'Median Bulk Probe'});
        title(sprintf('%% UMIs Retained=%5.2f', full(sum(counts.by_median_bulk_probe(:))/sum(umis_per_cell))*100));
    end
    
    subplot(2,4,4);
    if (isequal(map_method, 'max_probe_in_cell'))    
        semilogy([sort(reads_per_cell, 'descend') sort(sum(max_probe_reads,2), 'descend')]);
        axis tight;
        ylim([1 max(reads_per_cell)])
        xlabel('Cell (ranked');
        ylabel('Reads');    
        legend({'All Probes'; 'Max Probe'});        
        title(sprintf('%% Reads Retained=%5.2f', full(sum(max_probe_reads(:))/sum(reads_per_cell))*100));   
    else
        semilogy([sort(reads_per_cell, 'descend') sort(sum(bulk_median_probe_reads,2), 'descend')]);
        axis tight;
        ylim([1 max(reads_per_cell)])
        xlabel('Cell (ranked');
        ylabel('Reads');    
        legend({'All Probes'; 'Median Bulk Probe'});
        title(sprintf('%% Reads Retained=%5.2f', full(sum(bulk_median_probe_reads(:))/sum(reads_per_cell))*100));
    end
        
    subplot(2,4,[5, 6]);
    if (isequal(map_method, 'max_probe_in_cell'))    
        mask = counts.by_gene_max_probe>0;
        frac_express_max = arrayfun(@(i) sum(counts.which_probe_max(mask(:,i),i) == mode(counts.which_probe_max(mask(:,i),i)))/sum(mask(:,i)), [1:length(gene_name)]);    
        plot(sort(frac_express_max*100, 'descend'));
        axis tight;
        xlabel('Ranked Genes');
        ylabel('% Expressing Cells with Same Max Probe');
    else
        frac_express_bulk_median = mean(counts.by_median_bulk_probe>0,1);    
        plot(sort(frac_express_bulk_median*100, 'descend'));
        axis tight;
        xlabel('Ranked Genes');
        ylabel('% Cells Expressing Bulk Median');
    end
    
    subplot(2,4,[7,8]);    
    if (isequal(map_method, 'max_probe_in_cell'))    
        mask = counts.by_gene_total>0;
        histogram(counts.by_gene_max_probe(mask)./counts.by_gene_total(mask)*100)
        xlabel('% UMIs Retained Using Max Probe');
        ylabel('# of NonZero (Gene x Cell) Entries');
        title(sprintf('%% UMIs Retained=%5.2f', full(sum(counts.by_gene_max_probe(:))/sum(counts.by_gene_total(:))*100)));
    else
        mask = counts.by_gene_total>0;
        histogram(counts.by_median_bulk_probe(mask)./counts.by_gene_total(mask)*100);
        xlabel('% UMIs Retained Using Bulk Median');
        ylabel('# of NonZero (Gene x Cell) Entries');
        title(sprintf('%% UMIs Retained=%5.2f', full(sum(counts.by_median_bulk_probe(:))/sum(counts.by_gene_total(:))*100)));
    end

    suptitle(input_path);
    
    mkdir(outdir);
    print(sprintf('%s/QC.png', outdir), '-dpng', '-r600');
    close all;
    
    cb_path = sprintf('%s/barcodes.tsv', outdir);
    fid = fopen(cb_path, 'wt');
    fprintf(fid, '%s\n', CBs{:});
    fclose(fid);
    gzip(cb_path);
    delete(cb_path);
    
    gene_path = sprintf('%s/features.tsv', outdir);
    fid = fopen(gene_path, 'wt');
    for i = 1:length(gene_name)
        fprintf(fid, '%s\t%s\t%s\n', gene_name{i}, gene_name{i}, 'Gene Expression');
    end
    fclose(fid);
    gzip(gene_path);
    delete(gene_path);
    
    if (isequal(map_method, 'max_probe_in_cell'))    
        [r, c, v] = find(counts.by_gene_max_probe);
    else
        [r, c, v] = find(counts.by_median_bulk_probe);
    end
    mtx_path = sprintf('%s/matrix.mtx', outdir);
    fid = fopen(mtx_path, 'wt');    
    
    fprintf(fid, '%s\n', '%%MatrixMarket matrix coordinate integer general');
    fprintf(fid, '%s\n', '%metadata_json: {"format_version": 2, "software_version": "3.1.0"}');
    
    fprintf(fid, '%d %d %d\n', length(gene_name), length(CBs), length(v));
    for i = 1:length(v)
        fprintf(fid, '%d %d %d\n', c(i), r(i), v(i));
    end
    fclose(fid);
    gzip(mtx_path);
    delete(mtx_path);
    
    save(sprintf('%s/Transcriptome.mat', outdir), 'counts', 'dat', 'CBs', 'gene_name', 'input_path', 'map_method');
    
end