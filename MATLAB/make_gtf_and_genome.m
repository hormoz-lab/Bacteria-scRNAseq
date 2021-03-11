function make_gtf_and_genome(input_file, outdir, sample)
%make_gtf_and_genome creates GTF and FASTA files from an input table of
%probes
%
%   make_gtf_and_genome(INPUT_FILE, OUTDIR, SAMPLE) creates SAMPLE.FASTA and
%   SAMPLE.GTF in OUTDIR. INPUT_FILE contains 3 columns with each row 
%   corresponding to a probe. The 1st column is the index of a probe for a 
%   particular gene. The 2nd columns is the gene name. The 3rd column is the 
%   probe nucleotide sequence.
%
%   Author: Duluxan Sritharan. Hormoz Lab. Harvard Medical School.

    assert(exist(input_file)==2, 'Missing input file: %s', input_file);

    probe_file = table2cell(readtable(input_file));
    if (isnan(probe_file{end,1}))
        probe_file = probe_file(1:end-1,:);
    end
    [~, ~, which] = unique(probe_file(:,3), 'stable');    
    probe_file(ismember(which, find(accumarray(which,1)>1)),:) = [];
    
    [gene, ~, which] = unique(probe_file(:,2), 'stable');    
    
    mkdir(outdir);
    fid_gtf   = fopen(sprintf('%s/%s.gtf', outdir, sample),   'wt');
    fid_fasta = fopen(sprintf('%s/%s.FASTA', outdir, sample), 'wt');
    fprintf(fid_fasta, '>chr1\n');
    
    s = 0;
    
    for i = 1:length(gene)
        probes = find(which==i);
        assert(all(diff(probes)==1));
        for j = 1:length(probes)
            probe_seq = probe_file{probes(j),3};
            fprintf(fid_fasta, '%s', probe_seq);            
            probe_name = sprintf('%s_%d', gene{i}, j);
            
            fprintf(fid_gtf, '%s\t%s\t%s\t%d\t%d\t%c\t%c\t%d\t%s "%s"; %s "%s";\n', ...
                    'chr1', probe_name, 'exon', s+1, s+length(probe_seq), '.', '+', 0, 'gene_id', probe_name, 'transcript_id', probe_name);

            s = s+length(probe_seq);
        end
    end
    fprintf(fid_fasta, '\n');
    
    fclose(fid_fasta);
    fclose(fid_gtf);
    
    gtf = splitlines(fileread(sprintf('%s/%s.gtf', outdir, sample)));
    gtf = gtf(1:end-1);
    gtf = cellfun(@strsplit, gtf, 'un', false);
    gtf = vertcat(gtf{:});
    gene_s = cellfun(@str2num, gtf(:,4));
    gene_e = cellfun(@str2num, gtf(:,5));
    assert(isequal(gene_s(2:end),gene_e(1:end-1)+1));
    
    probe_name = gtf(:,2);
    assert(length(probe_name)==length(unique(probe_name)));
    gene_name = cellfun(@(x) strsplit(x, '_'), probe_name, 'un', false);
    gene_name = cellfun(@(x) horzcat(x{1:end-1}), gene_name, 'un', false);
    assert(isequal(gene_name, probe_file(:,2)));
    
    genome = splitlines(fileread(sprintf('%s/%s.FASTA', outdir, sample)));
    genome = horzcat(genome{2:end});
    probe_seqs = arrayfun(@(s,e) genome(s:e), gene_s, gene_e, 'un', false);
    assert(length(probe_seqs)==length(unique(probe_seqs)));
    
    assert(isequal(probe_seqs, probe_file(:,3)));
    
end