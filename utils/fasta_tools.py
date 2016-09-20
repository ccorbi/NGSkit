import os


def write_fasta_sequence(sequence_info, file_output, write_mode='a'):
    """Add sequences to a file, in Fasta Format.

    Parameters
    ----------
    sequence_info : lisrt, tuple
        sequence_info[0] == Header or id of the sequences, if do not contain > ,
    it will be added.
        sequence_info[1] == Sequence



    Returns
    -------
        let the file_handle open, this should improve I/O


    """
    # check input lenght
    assert(len(sequence_info)>=2)

    if isinstance(file_output, str):
        file_handle = open(file_output, write_mode)
    else:
        file_handle = file_output

    file_handle.write(">{0}\n{1}\n".format(*sequence_info))


    return


def to_fasta(grp_seq, output='file.fa', header=False):

    if isinstance(grp_seq, dict):
        dict_to_fasta(grp_seq)
    if isinstance(grp_seq, pd.DataFrame):
        df_to_fasta(grp_seq)
    if isinstance(grp_seq, str):
        write_seq(read1_id, read1_seq, suffix, l, out_dir)
    else:
        raise TypeError
