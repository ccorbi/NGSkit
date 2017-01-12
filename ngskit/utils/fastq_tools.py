"""FastQ tools
"""
def write_fastq_sequence(sequence_info, file_output, write_mode='a'):
    """Add sequences to a file, in Fastq Format.

    Parameters
    ----------

    sequence_data : array_like
        sequence_info[0] == Header or id of the sequences, if do not contain @ ,
    it will be added.
        sequence_info[1] == Sequence
        sequence_info[2] == Quality

    file_output: str, obj
        This function can recive both a file_handler or file name. In the former
        scenario it will create a file_handler, and in both cases it will let
        it open, to improve I/O.


    Returns
    -------

    file_handle : obj
        returns the file handler.

    Raises
    ------
    ValueError
        Sequence_data should contain three  items: header, Sequence and PhedQ

    Examples
    --------


        >>> write_fasta_sequence(['SEQ_1',
                                  'ATGATGATGA',
                                  '/66AA<<<6A'],
                                  'my_file.fasta')

        >>> write_fasta_sequence(['SEQ_1',
                                  'ATGATGATGA',
                                  '/66AA<<<6A'],
                                  open('my_file.fasta', 'a'))


    """
    # check input lenght
    if not len(sequence_info)>=2:
        raise ValueError("Sequence data must contain at least header, sequence and Quality")

    if '@' in sequence_info[0]:
        sequence_info[0] = sequence_info[0].replace('@','')

    if isinstance(file_output, str):
        file_handle = open(file_output, write_mode)
    else:
        file_handle = file_output

    file_handle.write(
        "@{0}\n{1}\n+\n{2}\n".format(*sequence_info))

    return file_handle


def avg_score(qual):
    """Return AVG score of a Fastq sequence

    Parametres
    ----------
    qual : str
        string contain the quality of the sequence encode in Phred

    Returns
    -------
    float
        Average quality of the sequences

    """
    # Translate the Quality to numbers
    avgQual = [ord(c)-33 for c in qual.rstrip("\n")]

    return avgQual
