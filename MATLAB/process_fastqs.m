function process_fastqs(input_fastq, outdir, tenX_version)
%process_fastqs splits Illumina sequencing runs of bacterial scRNAseq into
%a format recognized by CellRanger using both 10X and Probe UMIs
%
%   process_fastqs(INPUT_FASTQ, OUTDIR, TENX_VERSION) splits INPUT_FASTQ
%   into R1 and R2 reads appropriate for CellRanger and places the resulting 
%   FastQ files in OUTDIR/10X (10X UMIs in R1) and OUTDIR/Probe
%   (Probe UMIs in R1). Specify the TENX_VERSION used to 
%   prepare the library (2 or 3), to set the UMI length trimmed from the
%   INPUT_FASTQ.
%
%   Author: Duluxan Sritharan. Hormoz Lab. Harvard Medical School.

    tic;
    
    assert(exist(input_fastq, 'file')==2);
    
    assert(~(exist(outdir, 'dir')==7), 'Output path already exists');
    mkdir(outdir);
    mkdir(sprintf('%s/10X', outdir));
    mkdir(sprintf('%s/Probe', outdir));    
    
    fprintf('Processing: %s\n', input_fastq);
    
    Extender = seqrcomplement('CATAGTTTCTTCGAGCAA');
    L_ext = length(Extender);
    L_CB = 16;    
    if (tenX_version==2)
        L_UMI = 10;
    elseif (tenX_version==3)
        L_UMI = 12;
    else
        error('Unrecognized 10X version');
    end
    
    [input_fastq, ext] = maybe_unzip(input_fastq);
    info = fastqinfo(input_fastq);
    N_reads = double(info.NumberOfEntries);
    block_size = 1e6;
    N_blocks = ceil(N_reads/block_size);
    fprintf('Splitting %10d reads over %4d blocks\n', N_reads, N_blocks);
    
    for b = 1:N_blocks
        tic;        
        lo = (b-1)*block_size+1;
        hi = min(b*block_size, N_reads);        
        [L_orig{b}, masks{b}, probe_umi_end_pos{b}, sc{b}] = ...
            process_per_block(input_fastq, lo, hi, Extender, L_CB, L_UMI, L_ext, outdir);
        fprintf('Block: %4d/%4d. Reads: [%10d, %10d] in %6d seconds\n', b, N_blocks, lo, hi, toc);
    end   

    L_orig = vertcat(L_orig{:});
    masks = vertcat(masks{:});
    probe_umi_end_pos = vertcat(probe_umi_end_pos{:});
    sc = vertcat(sc{:});
    
    maybe_clear_unzipped(input_fastq, ext);
    
    fprintf('Output R2...\n');    
    
    output_fastq2 = sprintf('%s/10X/R2_001.fastq', outdir);
    info = fastqinfo(output_fastq2);
    assert(double(info.NumberOfEntries) ==  length(vertcat(masks.valid_seq)));
    gzip(output_fastq2);
    delete(output_fastq2);
    copyfile(sprintf('%s.gz', output_fastq2), strrep(sprintf('%s.gz', output_fastq2), '10X', 'Probe'));
    
    fprintf('Output 10X R1...\n');    
    
    output_fastq1 = sprintf('%s/10X/R1_001.fastq', outdir);
    info = fastqinfo(output_fastq1);
    assert(double(info.NumberOfEntries) == length(vertcat(masks.valid_seq)));    
    gzip(output_fastq1);
    delete(output_fastq1);
    
    fprintf('Output Probe R1...\n');    
    
    output_fastq1 = sprintf('%s/Probe/R1_001.fastq', outdir);
    info = fastqinfo(output_fastq1);
    assert(double(info.NumberOfEntries) == length(vertcat(masks.valid_seq)));
    gzip(output_fastq1);
    delete(output_fastq1);
    
    fprintf('Output provenance...\n');    
       
    metric_file = sprintf('%s/meta.mat', outdir);
    save(metric_file, 'input_fastq', 'L_orig', 'L_CB', 'L_UMI', 'L_ext', 'masks', 'sc', 'probe_umi_end_pos', '-v7.3', '-nocompression');
    
    fprintf('%f seconds\n', toc);
    
end

function [L_orig, masks, probe_umi_end_pos, sc] = ...
    process_per_block(input_fastq, lo, hi, Extender, L_CB, L_UMI, L_ext, outdir)
    
    fprintf('\tReading S...\n');
    [~, S, ~] = fastqread(input_fastq, 'Blockread', [lo, hi]);
    S=S';
    
    fprintf('\tTrimming by length...\n');    
    L_orig = cellfun(@length, S);
    masks.full_seq = L_orig > 110;    
    S = S(masks.full_seq);    
    masks.full_seq = find(masks.full_seq);
    
    fprintf('\tAligning...\n');    
    [sc, al] = nwalign(S, Extender, 'Alphabet', 'NT', 'glocal', true);
    masks.extender_found = (sc>=17);
    al = al(masks.extender_found);    
    S = S(masks.extender_found);  
    masks.extender_found = masks.full_seq(masks.extender_found);
    
    fprintf('\tTrimming by extender...\n');    
    probe_umi_end_pos = cellfun(@(x) sum(x(1,1:find(x(3,:) ~= '-', 1, 'first')-1)~='-'), al);
    clear al;
    % Should be a polyT region separating separating 10X (CB/UMI) and Probe UMI
    probe_mask = find(probe_umi_end_pos >= L_CB+2*L_UMI);
    S = S(probe_mask);    
    masks.valid_probe = masks.extender_found(probe_mask);
    
    % Keep the extender since the sequences are ~30 bps without it which
    % may affect alignment, but remove sequences consisting of just the
    % extender
    seq_mask = (L_orig(masks.valid_probe)>probe_umi_end_pos(probe_mask)+L_ext);
    S = S(seq_mask);
    masks.valid_seq = masks.valid_probe(seq_mask);
    
    fprintf('\tSplitting S...\n');
    
    S_CB  = cellfun(@(x)   x(1:L_CB), S, 'un', false);
    S_10X_UMI = cellfun(@(x) x(L_CB+1:L_CB+L_UMI), S, 'un', false);
    S_probe_UMI = cellfun(@(x,e) x(e-L_UMI+1:e), S, num2cell(probe_umi_end_pos(probe_mask(seq_mask))), 'un', false);
    S = cellfun(@(x,s) x(s+1:end), S, num2cell(probe_umi_end_pos(probe_mask(seq_mask))), 'un', false);
    
    fprintf('\tReading Q...\n');
    [~, ~, Q] = fastqread(input_fastq, 'Blockread', [lo, hi]);
    Q=Q(masks.valid_seq)';
    
    fprintf('\tSplitting Q...\n');
    
    Q_CB  = cellfun(@(x)   x(1:L_CB), Q, 'un', false);
    Q_10X_UMI = cellfun(@(x) x(L_CB+1:L_CB+L_UMI), Q, 'un', false);
    Q_probe_UMI = cellfun(@(x,e) x(e-L_UMI+1:e), Q, num2cell(probe_umi_end_pos(probe_mask(seq_mask))), 'un', false);
    Q = cellfun(@(x,s) x(s+1:end), Q, num2cell(probe_umi_end_pos(probe_mask(seq_mask))), 'un', false);
    
    fprintf('\tReading H...\n');
    [H, ~, ~] = fastqread(input_fastq, 'Blockread', [lo, hi]);
    H=H(masks.valid_seq)';
    
    fprintf('\tOutput R2...\n');    
    
    output_fastq2 = sprintf('%s/10X/R2_001.fastq', outdir);           
    fastqwrite(output_fastq2, H, S, Q);    
    
    fprintf('\tOutput 10X R1...\n');    
    
    output_fastq1 = sprintf('%s/10X/R1_001.fastq', outdir);        
    S = strcat(S_CB, S_10X_UMI);
    Q = strcat(Q_CB, Q_10X_UMI);
    fastqwrite(output_fastq1, H, S, Q);
        
    fprintf('\tOutput Probe R1...\n');    
    
    output_fastq1 = sprintf('%s/Probe/R1_001.fastq', outdir);
    S = strcat(S_CB, S_probe_UMI);
    Q = strcat(Q_CB, Q_probe_UMI);
    fastqwrite(output_fastq1, H, S, Q);    
    
end
