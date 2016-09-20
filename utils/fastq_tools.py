import os


def write_fastq_sequence(sequence_info, file_output, write_mode='a'):
    """Write a sequence to a file. Fastq format.

    Parameters
    ----------
    sequence_info : lisrt, tuple
        sequence_info[0] == Header or id of the sequences, if do not contain @ ,
    it will be added.
        sequence_info[1] == Sequence
        sequence_info[2] == Quality


    Returns
    -------
     it lets the file_handle open. this should improve I/O

    """
    # check input lenght
    assert(len(sequence_info)>=3)

    if '@' in sequence_info[0]:
        sequence_info[0] = sequence_info[0].replace('@','')

    if isinstance(file_output, str):
        file_handle = open(file_handle, write_mode)
    else:
        file_handle = file_output

    file_handle.write(
        "@{0}\n{1}\n+\n{2}\n".format(*sequence_info))

    return


def get_score():
    pass
    return


def avg_score(read1_qual):

    # Translate the Quality to numbers
    qual = [ord(c)-33 for c in quality.rstrip("\n")]

    return
