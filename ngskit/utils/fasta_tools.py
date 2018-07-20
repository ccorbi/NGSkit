"""Fasta Tools
"""
def write_fasta_sequence(sequence_data, file_output, write_mode='a'):
    """Add sequences to a file, in Fasta Format.


    Parameters
    ----------
    sequence_data : str
        Sequence to add to the fasta file. if only the sequence is provided,
    assume the header is not relevant and a random will be created, sequence base
    to avoid collisions

    sequence_data : array_like
        sequence_data[0] == Header or id of the sequences, if do not contain > ,
    it will be added.
        sequence_data[1] == Sequence

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
        Sequence_data should contain two items: header, Sequece
    Examples
    --------

        >>> write_fasta_sequence('ATGATGATGA','my_file.fasta')

        >>> write_fasta_sequence('ATGATGATGA',open('my_file.fasta', 'a'))

        >>> write_fasta_sequence(['SEQ_1', 'ATGATGATGA'],'my_file.fasta')


    """
    # Check the input sequence
    if isinstance(sequence_data, str):
        # create a  Header using 100 first sequence caracters.
        header = sequence_data.strip('\n').strip()[:100]
        sequence_data = [header,
                        sequence_data.strip('\n').strip()]

    if not len(sequence_data)>=2:
        raise ValueError("Sequence data must contain at least header and sequence")

    # check if a file handelr has been provided
    if isinstance(file_output, str):
        file_handle = open(file_output, write_mode)
    else:
        file_handle = file_output

    # write the sequence
    file_handle.write(">{0}\n{1}\n".format(*sequence_data))

    return file_handle


def to_fasta(grp_seq, output, header=False):
    """Transform a batch of sequnces to a fasta format file.

    Parameters
    ----------

    grp_seq : array_like
        Iterable object with sequneces


    """

    if header == False:
        for sequence in grp_seq:
            output = write_fasta_sequence(sequence, output, write_mode='a')
        output.close()



        
        
def read_fasta(filn,proc_header=None):
    
    default_info = dict()
    header = None
    seq = list()
    
    with open(filn, 'r') as fn:
        for line in fn:
            if line[0] == '>':
                    if proc_header:
                        header = proc_header(header)
                        
                    if header:
                        default_info[header] = ''.join(seq)
                    header = line[1:].strip()
                    seq = list()
                                            
            else:
                seq.append(line.strip())
    
    default_info[header] = ''.join(seq)
    return default_info

def read_to(filename, to='dict', proc_header=None):
    
    default = read_fasta(filename, proc_header)
    
    if to == 'dict':
        
        return default
    
    elif to == 'list':
        return list(default.values())
    
    elif to == 'pandas':
        data = zip(default_info.keys(),default_info.values())
        return pd.DataFrame(data, columns = ['ID','Seq'])
        