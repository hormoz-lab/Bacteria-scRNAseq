function [counts, barcodes, genes] = load_10x_mtx(folder)

    assert(exist(folder, 'dir')==7, sprintf('Invalid path to 10X Matrix folder: %s\n', folder));
    
    fprintf('Loading 10X Matrix folder: %s\n', folder);
    
    if (exist(sprintf('%s/barcodes.tsv.gz', folder), 'file')==2)
        gunzip(sprintf('%s/barcodes.tsv.gz', folder));
        barcodes = splitlines(fileread(sprintf('%s/barcodes.tsv', folder)));
        delete(sprintf('%s/barcodes.tsv', folder));
    elseif (exist(sprintf('%s/barcodes.tsv', folder), 'file')==2)
        barcodes = splitlines(fileread(sprintf('%s/barcodes.tsv', folder)));
    end
    barcodes = barcodes(1:end-1);    
    
    if (exist(sprintf('%s/features.tsv.gz', folder), 'file')==2)
        gunzip(sprintf('%s/features.tsv.gz', folder));
        genes = splitlines(fileread(sprintf('%s/features.tsv', folder)));
        delete(sprintf('%s/features.tsv', folder));
    elseif (exist(sprintf('%s/features.tsv', folder), 'file')==2)
        genes = splitlines(fileread(sprintf('%s/features.tsv', folder)));
    end
    genes = cellfun(@(x) strsplit(x, '\t'), genes(1:end-1), 'un', false);
    genes = vertcat(genes{:});
    genes = genes(:,1);
    
    if (exist(sprintf('%s/matrix.mtx.gz', folder), 'file')==2)
        gunzip(sprintf('%s/matrix.mtx.gz', folder));
        counts = splitlines(fileread(sprintf('%s/matrix.mtx', folder)));
        delete(sprintf('%s/matrix.mtx', folder));
    elseif (exist(sprintf('%s/matrix.mtx', folder), 'file')==2)
        counts = splitlines(fileread(sprintf('%s/matrix.mtx', folder)));
    end
    counts = cellfun(@str2num, counts(3:end-1), 'un', false);
    counts = vertcat(counts{:});    
    meta = counts(1,:);
    counts = counts(2:end,:);
    assert(max(counts(:,1))<=meta(1));
    assert(max(counts(:,2))<=meta(2));
    assert(length(counts)  ==meta(3));
    counts = sparse(counts(:,2), counts(:,1), counts(:,3), meta(2), meta(1));
    
    assert(size(counts,1) == size(barcodes,1), 'Number of barcodes and rows in count matrix do not match');
    assert(size(counts,2) == size(genes,1), 'Number of genes and columns in count matrix do not match');
    
end

