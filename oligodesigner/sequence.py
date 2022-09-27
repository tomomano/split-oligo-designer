import mygene
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd

def fetch_mygene(genename, spiecies):
    genename = genename.upper()
    mg = mygene.MyGeneInfo()
    geneinfo = mg.querymany(genename, scopes='symbol', fields='ensembl.gene', species=spiecies,
                            verbose=False)
    return geneinfo


def run_blast(input_file, output_file, database):
    """
    run blast against transcriptome database. negative_seqid option does not work for an unknown reason.
    """
    blastn_cline = NcbiblastnCommandline(
        task='blastn-short',
        query=input_file,
        db=database,
        # num_threads = 4,
        strand='minus',
        outfmt='10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        out=output_file
    )
    stdout, stderr = blastn_cline()


def exclude_self(df_blast, geneinfo):
    """
    omit the gene of interest from the blast result
    by referring the gene information from mygene.
    """
    mg = mygene.MyGeneInfo()
    transcript_series = df_blast['sseqid'].str.split('.').str[0].drop_duplicates(keep='first')
    df_transcript_id = pd.DataFrame(mg.querymany(transcript_series, scopes='ensembl.transcript', fields='id', verbose = False))
    gene_id = geneinfo[0]['_id']
    df_selfrmv = df_blast[~df_blast['sseqid'].str.split('.').str[0].isin(df_transcript_id[df_transcript_id['_id']==gene_id]['query'])]
    return df_selfrmv.reset_index(drop=True)


def get_off_targeting_oligo(minimum_offtarget_gap, df_blast):
    """

    """
    oligo_to_be_removed = []
    df_nonspecific_oligo = df_blast[df_blast['sseqid'].duplicated(keep=False)]
    # check each oligo whether the oligo potentially causes a off-target signal
    for nonspecific_oligo in list(df_nonspecific_oligo['qseqid'].unique()):
        nonspecific_gene_set = (df_nonspecific_oligo[df_nonspecific_oligo['qseqid'] == nonspecific_oligo])['sseqid']
        # check each nonspecific binding site.
        for nonspecific_gene in list(nonspecific_gene_set.unique()):
            # calculate distance to the nearest off-target binding site.
            sstart_send = df_nonspecific_oligo[df_nonspecific_oligo['sseqid'] == nonspecific_gene]
            sstart_send = (sstart_send.set_index('qseqid', drop=True))[['sstart', 'send']]
            sstart_send = sstart_send.sort_values(by=['sstart'], ascending=True)
            dist_before = (sstart_send['send'] - sstart_send['sstart'].shift(periods=1)).rename('dist_before')
            sstart_send = pd.concat([sstart_send, dist_before, (dist_before.shift(periods=-1)).rename('dist_after')],
                                    axis=1)
            # see whether the binding site is within the intolerable range.
            leaky_set = sstart_send[(sstart_send['dist_before'].abs() < minimum_offtarget_gap) | (
                        sstart_send['dist_after'].abs() < minimum_offtarget_gap)]

            # if oligo is causing off-target signal, remove.
            if nonspecific_oligo in list(leaky_set.index):
                df_nonspecific_oligo = df_nonspecific_oligo[~(df_nonspecific_oligo['qseqid'] == nonspecific_oligo)]
                # update the oligo list
                oligo_to_be_removed.append(nonspecific_oligo)
    return oligo_to_be_removed


def remove_oligo(oligominer_df, oligo_to_be_removed):
    """
    """
    selected_oligo_df = oligominer_df[~(oligominer_df['oligo_id'].isin(oligo_to_be_removed))]
    selected_oligo_df = selected_oligo_df.reset_index(drop=True)
    return selected_oligo_df


def calc_oligo_interval(oligo_df, oligo_number):
    """
    """
    dist_before = (oligo_df['start'] - oligo_df['end'].shift(periods=1)).rename('interval_before')
    dist_after = (oligo_df['start'].shift(periods=-1) - oligo_df['end']).rename('interval_after')
    # To select the probe sets, scan the sum of gaps between each pair in 24 sets, then find minimal point.
    mean_interval = dist_before.rolling(oligo_number, min_periods=oligo_number).mean()
    mean_interval = mean_interval.rename('mean_interval')
    interval_df = pd.concat([dist_after, mean_interval], axis=1)
    return interval_df


def add_hcr_seq_v3(selected_oligo_df, hcr_initiator_set):
    """
    add HCR sequences to both sides of split probes.
    """
    splitseq_ll = hcr_initiator_set['splitseq_ll']
    splitseq_rr = hcr_initiator_set['splitseq_rr']
    splitseq_lr = hcr_initiator_set['splitseq_lr']
    splitseq_rl = hcr_initiator_set['splitseq_rl']
    anchorseq_ll = hcr_initiator_set['anchorseq_ll']
    anchorseq_rr = hcr_initiator_set['anchorseq_rr']
    anchorseq_lr = hcr_initiator_set['anchorseq_lr']
    anchorseq_rl = hcr_initiator_set['anchorseq_rl']

    selected_oligo_df = selected_oligo_df.sort_values(by='start')

    even_oligo = selected_oligo_df[selected_oligo_df.index % 2 == 0]
    even_hcr_oligo = (splitseq_ll + anchorseq_ll + even_oligo['seq'] + anchorseq_lr + splitseq_lr).rename('hcr_seq')
    even_oligo = pd.concat([even_oligo, even_hcr_oligo], axis=1)

    odd_oligo = selected_oligo_df[selected_oligo_df.index % 2 == 1]
    odd_hcr_oligo = (splitseq_rl + anchorseq_rl + odd_oligo['seq'] + anchorseq_rr + splitseq_rr).rename('hcr_seq')
    odd_oligo = pd.concat([odd_oligo, odd_hcr_oligo], axis=1)

    hcr_oligo_set = (pd.concat([even_oligo, odd_oligo], axis=0)).sort_values(by='start')
    return hcr_oligo_set