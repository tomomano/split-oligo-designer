from oligodesigner import blockParse_py3
from oligodesigner import io
from oligodesigner import parse
from oligodesigner import sequence
import pandas as pd


def generator(mFISH3D_param, oligominer_param, num_threads=1):
    """
    Arguments:
        mFISH3D_param (dict)
            'fasta': fasta file of the template sequence
            'database': database for blast
            'minimum_offtarget_gap': If the gap between two non-specific binding is more than minimum_offtarget_gap,
                the pair is not regarded to cause a off-target signal. Recommended value: 100.
            'hcr_seqs': sequences of HCR fragments. The example design can be found in oligodesigner/parameters.py
            'self_remove': Recommended: True. set self_remove True if you want to remove the template sequence from
                off-target analysis. Otherwise, the gene of your interest could be regarded as an off-target product.
                Turn this to False when your template is not found in database (e.g. GFP).
    """
    fasta = mFISH3D_param['fasta']
    database = mFISH3D_param['database']
    minimum_offtarget_gap = mFISH3D_param['minimum_offtarget_gap']
    hcr_initiator_set = mFISH3D_param['hcr_seqs']
    self_remove = mFISH3D_param['self_remove']

    # Make path to temp files and result file
    (oligominer_fastq, oligominer_fasta, oligo_sets_path, param_path) = io.gen_file_path(fasta)

    # Save parameters to csv file.
    io.record_param(mFISH3D_param, oligominer_param, param_path)

    # Run OligoMiner
    blockParse_py3.runSequenceCrawler(fasta, *parse.oligominer_param_parser(oligominer_param))

    # read oligominer fastq as dataframe
    oligominer_df = io.dataframe_from_oligominer(oligominer_fastq)

    # converting fastq to fasta. writing fasta is a must for blast.
    io.convert_fastq2fasta(oligominer_fastq, oligominer_fasta)
    io.add_id(oligominer_fasta)

    if self_remove:
        seqid = sequence.refseqid_from_fasta(fasta)
        # seqid = sequence.check_refseqid_exist(fasta, return_entrezid=False)
        # if seqid is None:
        #     seqid = sequence.get_homology_in_database(fasta, database, num_threads=num_threads)
        # else:
        #     _ = sequence.get_homology_in_database(fasta, database, num_threads=num_threads)

        # run blast
        df_blast = sequence.run_blast_df(oligominer_fasta, database, task='blastn-short', strand='plus', num_threads=num_threads)

        # Exclude sequences that are homologous to the gene of a interest from the result
        df_blast = sequence.exclude_self(df_blast, seqid)

        # Remove combination if sstart and send are close enough and if they are rective combination.
        oligo_to_be_removed = sequence.get_off_targeting_oligo(minimum_offtarget_gap, df_blast)
        print("These sequences were removed because BLAST detected potential off-target binding", oligo_to_be_removed)

        oligo_df = sequence.remove_oligo(oligominer_df, oligo_to_be_removed)

    else:
        df_blast = sequence.run_blast_df(oligominer_fasta, database, task='blastn-short', strand='plus', num_threads=num_threads)
        oligo_to_be_removed = sequence.get_off_targeting_oligo(minimum_offtarget_gap, df_blast)
        oligo_df = sequence.remove_oligo(oligominer_df, oligo_to_be_removed)

    oligo_df = sequence.reverse_complement_dataframe(oligo_df)

    interval = sequence.calc_oligo_interval(oligo_df)
    oligo_df = pd.concat([oligo_df, interval], axis=1)

    # Make oligos with HCR sequences.
    res = sequence.add_hcr_seq_v3(oligo_df, hcr_initiator_set)
    res.to_csv(oligo_sets_path, sep=',', index=False, header=True)

    return res

