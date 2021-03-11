function out = parse_molecule_info(h5_file, CBs, N_genes)

    assert(exist(h5_file, 'file')==2, 'Invalid path to molecule_info file')

    barcode_list      = h5read(h5_file, '/barcodes');
    [~, barcode_mask] = ismember(barcode_list, CBs);
    clear barcode_list

    barcode_id        = h5read(h5_file, '/barcode_idx')+1;
    barcode_mask      = barcode_mask(barcode_id);
    barcode_id        = nonzeros(barcode_mask);
    barcode_mask      = logical(barcode_mask);

    reads             = h5read(h5_file, '/count');
    reads             = reads(barcode_mask);

    gene_id           = h5read(h5_file, '/feature_idx')+1;
    gene_id           = gene_id(barcode_mask);

    umis              = h5read(h5_file, '/umi');
    [umis, ~, umi_id] = unique(umis(barcode_mask), 'rows');
    
    out.barcode_id    = barcode_id;
    out.read          = reads;
    out.gene_id       = gene_id;
    out.umis          = umis;
    out.umi_id        = umi_id;
    
    assert(size(unique([out.barcode_id out.umi_id], 'rows'),1) == length(out.gene_id));
   
    out.umi_counts  = accumarray([out.barcode_id, out.gene_id], 1, [length(CBs), N_genes], [], [], true);
    out.read_counts = accumarray([out.barcode_id, out.gene_id], double(out.read), [length(CBs), N_genes], [], [], true);
    
end