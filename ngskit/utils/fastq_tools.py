"""Common set of tools for FastQ files
"""
import gzip

def read_fastq(demultiplexed_fastq, qthreshold=0):
    """Read a fastq file, return a list of tuples seq:average_quality:quality, 
    also filter out sequence with an average quality below  qthreshold"""
    seqs = list()

    # Q&D hack
    if '.gz' in demultiplexed_fastq:
        o = gzip.open
    else:
        o = open

    with o(demultiplexed_fastq, 'r') as read1:
        for read1_id in read1:
            # Read 4 by 4
            # ID lane info, seq info etc
            # Read seq and Quality info
            read1_seq, read1_strand, read1_qual = [next(read1) for _ in range(3)]
            #Translate the Quality to a list of Integers
            try:
                qual = [ord(c)-33 for c in read1_qual.strip().decode()]
            except AttributeError:
                qual = [ord(c)-33 for c in read1_qual.strip()]

            # Control
            try:
                avg_quality = sum(qual)/float(len(qual))
            except ZeroDivisionError:
                print('Sequence with no lenght or no score')
                print(read1_id,read1_seq,read1_qual)
                raise ValueError

            if len(read1_seq.strip()) == len(qual) and avg_quality >= qthreshold:
                try:
                    seqs.append((read1_seq.strip().decode(), avg_quality, qual))
                except AttributeError:
                    seqs.append((read1_seq.strip(), avg_quality, qual))
            else:
                # Stats
                pass

    return seqs

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
