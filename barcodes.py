from __future__ import print_function
import os
import logging

logger = logging.getLogger(__name__)

class Barcode(object):
    """docstring for Barcode Class

    """
    def __init__(self, name, arg):
        """init barcodes. from a file.


        Parameters
        ----------
        name : str
            Name of the sample, or poll of sequences
        arg : list
            list of str, with b1, c1, c2, b2 and len

        Attributes
        ----------
        id : str
          sample name, or id
        b1_seq : str
          Seq first barcode
        c1_seq : str
          Seq Constant region 1
        c2_seq : str
          Seq Constant region 2
        b2_seq : str
          Seq barcode 2

        trgt_len : int
          Length, target sequence

        b1_len : int
          Length seq first barcode
        c1_len : int
          Length seq Constant region 1
        c2_len : int
          Length seq Constant region 2
        b2_len : int
          Length seq barcode 2


        """

        self.id = name
        self.b1_seq = arg[0].strip()
        self.c1_seq = arg[1].strip()
        self.c2_seq = arg[2].strip()
        self.b2_seq = arg[3].strip()
        self.trgt_len = int(arg[4].strip())
        self._calc_lens()

    def _calc_lens(self):
        self.b1_len = len(self.b1_seq)
        self.c1_len = len(self.c1_seq)
        self.b2_len = len(self.b2_seq)
        self.c2_len = len(self.c2_seq)

        return


def read(barcode_file):
    '''Method to Read barcode file, all the seq must be in 5' 3' sense.

    Parameters
    ----------

    barcode_file :  str
      barcode_file: The path or name of the file that contains barcode info

    Returns
    -------

    dict
      dictionary with barcode objects, name is the sample_id
      for a list of sequences.


    '''
    assert os.path.isfile(barcode_file)

    barcodes = list()

    with open(barcode_file, 'r') as input_file:
        for line in input_file:
            line = line.strip()
            # skip comments
            if not line.startswith('#'):
                data = line.split()
                # Sample ID
                if len(data) >= 6:
                    name = data[0].strip()
                    # Forward-Barcode, Barcode-2-Reversed, Konstant_region1, Konstant_region2-Reversed
                    barcodes.append(Barcode(name, [data[1].strip(), data[2].strip(),
                                                   data[3].strip(), data[4].strip(),
                                                   data[5].strip()]))

                    logger.info('BARCODE {} > b1:{} c1:{} c2:{} b1:{} target:{}'.format(name,
                                                      data[1].strip(),
                                                      data[2].strip(),
                                                      data[3].strip(),
                                                      data[4].strip(),
                                                      data[5].strip()))
                elif len(data) == 4:
                    barcodes.append(Barcode(name, [data[1].strip(), data[2].strip(),
                                                   '-', '-', data[3].strip()]))

                    logger.info('BARCODE {} > b1:{} c1:{} target:{}'.format(name,
                                                      data[1].strip(),
                                                      data[2].strip(),
                                                      data[3].strip()))
                else:
                    print('''Barcode should contain at least:\n
                               Sample ID, Forward-Barcode, Constant_region 1\n
                               and ideally:
                               Sample ID, Forward-Barcode, Constant_region 1,
                               Barcode-2, Constant_region2''')
                    print(data)

    return barcodes
