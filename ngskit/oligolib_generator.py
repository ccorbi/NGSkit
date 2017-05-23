# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 11:19:51 2015
Small script and library of functions to help to create the in house
oligochips.

@author: ccorbi


Example
-------

    > from  oligolib_generator import *
    > disrd_CONSTANT_F = 'CAGCCTCTTCATCTGGC'
    > disrd_CONSTANT_R = 'GGTGGAGGATCCGGAG'
    > peptides = {'RTLSQLYCSYGPLT': 'MUS81',
                   'NPFREKKFFCTIL': 'GNG4',
                   'TEGPDSD': 'SFN',
                   'PSSLAYSLKKH': 'INCEP',
                   'ETFSDLW': 'MDM2'}

    > for seq,_id in peptides.items():
    ...     lib = peptide_library(template=seq,
                                  CONSTANT_F = 'CAGCCTCTTCATCTGGC',
                                  CONSTANT_R = 'GGTGGAGGATCCGGAG')
    ...     lib.generate_single_variants()
    ...     lib.write('./fasta/{}_disorderome_CR.fasta'.format(_id), extend=87)


Notes
-----

    Tested on python 2.7. 3.6

    To Do ::
        Finish stand alone, configration and arguments passing
        Improve Testing by using pytests, and logging




"""
from __future__ import print_function
import sys
import random
import itertools
import json
import Levenshtein as ltein
import logging
import argparse
import time

from ngskit.utils.dna import (translate2aa, translate2na, clean_restriction_sites)


class peptide_library(object):
    """Main Class to Generate library.


    """

    # compatible class arguments


    def __init__(self, template=None, **kwargs):
        """Init.

        Parameters
        ----------

            seq : str, optional
                Template sequence to build the library
                Use None when you want to load variables from a file.
                If a template is not provided, all the generation options
                are disabled, thought.

        Kwargs:
        -------

            lib_name : str, optional
                Name of the library
            CONSTANT_F : str
                Forward constant region
                first constant region 5'
            CONSTANT_R : str
                Reverse constant region
                second constant region 3'

            include_template : bool
                The original template will be introduced in the library,
                (default True)

            restriction_enzymes : dict or list
                dictionary or list with the restriction sites
                (by default, dict_restriction = ['GAATTC','CCCGGG']
                or ecorI, XmaI, salI

            species : str
                species of the lirbary. This determine the codon frequency usage,
                (by default ('species':'human') accept 'E.coli')

        See Also
        --------
            load_designs
            write
            change_template



        """
        self._class_arguments = ['include_template',
                                 'lib_name',
                                 'CONSTANT_F',
                                 'CONSTANT_R',
                                 'lib_limit',
                                 'restriction_enzymes',
                                 'codon_usage_species']

        # check the passed arguments
        for k in kwargs:
            if k not in self._class_arguments:
                raise ValueError("{} got an unexpected keyword argument '{}' ".format(self.__class__.__name__, k))


        # init template sequence
        if isinstance(template, str):
            self.template = Template(template)

            # Include template in the library
            self.include_template = kwargs.get('include_template', True)

            # Optional name of the lib, by default use the sequences
            self.lib_name = kwargs.get('lib_name', self.template.seq)

        else:
            if template is None:
                self.template = False
                # Include template in the library
                self.include_template = kwargs.get('include_template', False)
                self.lib_name = kwargs.get('lib_name', 'No_Template')
            else:
                raise Exception

        # Set up Constant regions
        self.CONSTANT_F = kwargs.get('CONSTANT_F', '')
        self.CONSTANT_R = kwargs.get('CONSTANT_R', '')

        # limit for the library
        self.lib_limit = kwargs.get('lib_limit', None)

        # SET UP Restriction enzymes
        self.restriction_enzymes = dict()
        self.add_restriction_enzyme(kwargs.get('restriction_enzymes', dict()))
                                                    #{'ecorI': 'GAATTC',
                                                    #'XmaI': 'CCCGGG',
                                                    #'salI': 'GTCGAC'}
        # Setup CODON usage and validate
        self.codon_usage_species = kwargs.get('codon_usage_species', 'human') # human E.coli

        if self.codon_usage_species not in ['human', 'E.coli']:
            raise ValueError("{} is not supported codon usage ".format(self.codon_usage_species))


        # Init Internal Variables

        self._aalibrary = dict()

        if self.include_template:

            index_id = len(self._aalibrary)
            self._aalibrary[template] = self.lib_name + '_OO_' + str(index_id)

        return

    def __call__(self):
        self.info()

    def __len__(self):

        return len(self._aalibrary)

    def info(self):

        for arg in self._class_arguments:
            print("{} :\t{}".format(arg, self.__getattribute__(arg)))

        n = len(self._aalibrary)
        print("{} :\t{}".format('Seq in the lib', self.__getattribute__(n)))

        return

    def load_designs(self, seq_file, **kwargs):
        """ Load a file with the seq of designs for the template. No headers
        By default overwrite identical previous designs.


        Parameters
        ----------
            seq_file : str
                Path and name of the sequence file

        Kwargs:
        -------
            sep : str, optional
                separation of columns (default None)
            column : int, optional
                column were the sequences are (default 0)
            suffix : srt, optional
                Suffix of the variants (default _MD_)

        Raises:
        ------
            ValueError
                If lenght of the sequence in the file is different
                from the template

        Returns
        -------

        """

        # Set up suffix for seq id purpuses
        # _MD_ for default
        # Manual Desings
        suffix = kwargs.get('suffix', '_MD_')
        # Load file cfg ToDo use pandas
        sep = kwargs.get('sep', None)
        column = kwargs.get('column', 0)

        # Open file
        with open(seq_file, 'r') as input_file:
            for line in input_file:
                # if the file is empty or is a comment, skip
                if len(line.strip()) == 0 or line.startswith('#'):
                    continue

                if sep != False:
                    try:
                        if sep == '':
                            seq = line.split()[column]
                        else:
                            seq = line.split(sep)[column]
                    except:
                        print('ERROR!: Reading File {}  on:'.format(seq_file))
                        print(line)
                        raise

                else:
                    seq = line.strip()


                if not seq in self._aalibrary:
                    # check library limit
                    if self._lib_notfull() is False:
                        self._aalibrary[seq] = self._name_seq(suffix=suffix)

        return

    def write(self, file_name='', add_stop_codon=True, extend=False, **kwargs):
        """Translate and Write the library in fasta format.

        Parameters
        ----------
            file_name : str
                Ouput file  name, by default library name.fasta
            stop_end : bool
                 Add a double stop codon after the design, by default True
            extend : int
                Add stop codons to extend the final oligo sequences
                to  the given legnth, the constant regions  are included by default No.


        Returns
        -------


        """
        if file_name == '':
            file_name = self.lib_name + '.fasta'

        output = open(file_name, 'w')

        CONSTANT_F = kwargs.get('CONSTANT_F', self.CONSTANT_F)
        CONSTANT_R = kwargs.get('CONSTANT_R', self.CONSTANT_R)

        # This can be improved
        SPACER = 'TGATAATAGTAATAGTGATAGTGATAATGATAATGA'

        for seq, seq_id in self._aalibrary.items():

            seq_dna = translate2na(seq, species=self.codon_usage_species)
            # check for the restriction sites
            seq_dna_checked = clean_restriction_sites(seq_dna, self.restriction_enzymes)

            # add double stop codons at the end
            if add_stop_codon:

                stop = 'TGATAA'
            else:
                stop = ''

            # Extend  Oligo to reach a concrete length

            if extend:

                # Determine the legnth of the spacer sequence
                base_oligo = len(seq_dna_checked) + len(CONSTANT_F) + len(CONSTANT_R) + len(stop)
                spacer_len = extend - base_oligo
                if spacer_len < 0:
                    raise ValueError(
                        'ERROR,  refill len is {} and the raw oligo is {}'.format(extend, base_oligo))

                refill_seq = SPACER[:spacer_len]

            else:

                refill_seq = ''
            # save in fasta format
            print('>{}_{}'.format(seq_id, seq), file=output)
            print(CONSTANT_F + seq_dna_checked.lower() +
                  stop + refill_seq + CONSTANT_R, file=output)

        output.close()

        return

    def change_template(self, seq):
        '''Change or add a the template '''

        self.template = Template(seq)

        return

    def add_restriction_enzyme(self, enzymes):
        """Add restriction enzyme definition.

        Restriction enzymes definitions are restrictions used when the AA library is converted
        to NA. see write fucntion.

        Parameters
        ----------
        enzymes : dict, list
            Restriction enzyme definition, dictionary  `name: target_sequence`, or if is a list
            provide only the `target_sequences` and a number will be provided as a name

        Returns
        -------

        """
        if isinstance(enzymes, dict):
            for name, enzym in enzymes.items():
                self.restriction_enzymes[name] = enzym
        if isinstance(enzymes, list):
            idx = len(self.restriction_enzymes)
            for name, enzym in enumerate(enzymes):
                self.restriction_enzymes[name+idx] = enzym

        return


    def random_purge(self, nsize):
        """Remove randon N items from the library.

        Parameters
        ----------
        nsize : int
            Number of random items to remove


        """
        for _ in range(nsize):
            self._aalibrary.pop(random.choice(list(self._aalibrary)))

        return

    def generate_single_variants(self, **kwargs):
        """From a sequence (string), return a dictionary of single muntatns.

        Parameters
        ----------
            :param positions : List of position to be modified
            :rtype: list

            :param inrange: [start,end] coordinates inside the seq by default full length

            :param bias: ['A','W'] aminoacids excluded or restrict
            :rtype list

            :param bias_type: how to handle the bias AA,  exclude or restrict , by default exclude
            :rtype str

        Returns
        -------

        Raises
        ------


        """
        if self.template is False:
            raise Exception('A template must be defined to use this method')

        # t = self.template
        # if postions is not given use the entire lenght of the template
        positions = kwargs.get('positions', range(self.template.len_))

        # Add suffix to this peptides, by default _SV_ SingleVariant
        suffix = kwargs.get('suffix', '_SV_')

        for eachposition in positions:
            for itera in self.template.get_all_mutants_in(eachposition, **kwargs):

                a = ''.join(itera)
                # if the peptide exist in the library ignore
                if not a in self._aalibrary:
                    # check library limit
                    if self._lib_notfull() is False:
                        self._aalibrary[a] = self._name_seq(suffix=suffix)

        return

    def generate_inframe_variants(self, frame_size=2, **kwargs):
        """From a sequence (string), return a dictionary of  muntats in frame

        Parameters
        ----------

            template_seq : str
                Template sequence (string)
            mutants_library : dict
                dictionary with the library
            template_name : str
                id name for the sequences (string)
            frame_size : int
                not working. (default 2)

            inrange : list, tuple
                [start,end] coordinates inside the seq by default full length
            bias : list
                list wiht ['A','W'] aminoacids excluded

        Returns
        -------


        """
        if not self.template:
            raise Exception('A template must be defined to use this method')

        t = self.template

        suffix = kwargs.get('suffix', '_FM_')

        for itera in t.frame_mutations(**kwargs):

            a = ''.join(itera)

            if not a in self._aalibrary:
                # check library limit
                if self._lib_notfull() is False:
                    self._aalibrary[a] = self._name_seq(suffix=suffix)


        return

    def generate_random_variants(self, how_many=1000, mutant_kind=[2], **kwargs):
        """From a sequence (string), return a dictionary of random muntats.

        Parameters
        ----------


        lib_limit : int, (default 1000)
            maximum size of the lib. mutamnts adde to the lib is current size - limit
        mutant_kind : list
            Kind of mutants, double, triple, etc (list, by default double mutants [2])
        inrange: list, tuple
            [start,end] coordinates inside the seq by default full length
        bias : list, tuple
            ['A','W'] aminoacids excluded

        Returns
        -------



        """
        # mutants = {2:'double',3:'triple'}

        # generate double, triple, etc   random mutants

        if not self.template:
            raise Exception('A template must be defined to use this method')

        t = self.template
        # Random Variant _RA_
        suffix = kwargs.get('suffix', '_RA_')

        counter = 0
        while counter < how_many:
            for i in mutant_kind:
                # Soft random, at least 50% of the wildtype
                random_mutation = t.soft_randomization(num_mutations=i, **kwargs)

                if not random_mutation in self._aalibrary:
                    # check library limit
                    if self._lib_notfull() is False:
                        self._aalibrary[random_mutation] = self._name_seq(suffix=suffix)
                        counter += 1


        return

    def single_json(self, jsonpath):

        json_data = open(jsonpath).read()
        cfg_data = json.loads(json_data)

        for position, actions in cfg_data.items():

            self.generate_single_variants(positions=[int(position)], **actions)

        return

    def permutations_json(self, jsonpath, **kwargs):

        json_data = open(jsonpath).read()
        cfg_data = json.loads(json_data)
        permutation = []

        # Permutation
        suffix = kwargs.get('suffix', '_PR_')

        for eachposition in range(self.template.len_):
            permutation.append(cfg_data.get(str(eachposition), self.template.seq[eachposition]))

        for i in itertools.product(*permutation):
            variant = ''.join(list(i))

            if not variant in self._aalibrary:
                # check library limit
                if self._lib_notfull() is False:
                    self._aalibrary[variant] = self._name_seq(suffix=suffix)

        return


    def _name_seq(self, suffix):
        # Get index
        index_id = len(self._aalibrary)
        # I use _ as split of items in the fasta header
        if '_' in self.lib_name:
            lib_name = self.lib_name.replace('_', '-')
        else:
            lib_name = self.lib_name
    # Set up suffix for seq id purpuses

        return lib_name + suffix + str(index_id)


    def _lib_notfull(self):
        # check library limit
        # _library_ limit_satuts
        if self.lib_limit:
            if len(self._aalibrary) >= self.lib_limit:
                # convert this to a raise
                print('WARNING!: Library limit reached -{}-'.format(self.lib_limit))
                return True
            else:
                return False
        else:
            return False

###########


class Template(object):

    def __init__(self, seq):
        """Template Class

        Paramenters
        -----------

            seq : str seq:
                Aminoacid sequence

        """
        # Sequence
        self.seq = seq
        # len of the seq
        self.len_ = len(seq)
        # One letter code Amino acid list
        self._aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

        return

    def __str__(self):
        return self.seq

    def mod_template(self, position, aa):
        """introduce a mutation in the seq in a determinate postion.
        Can accept a sinple of multimple subsititutions

        Paramenters
        -----------

        postion : int, list
            Position or list of positions to change. from: 0->x.

        aa : str
            amminoacid one letter code, or list of aaminoacid.
                       must be syncronized  with the positions


        Returns
        -------
            str, list
            sequence with the modifications


        """
        variant = list(self.seq)

        if isinstance(position, list):
            if isinstance(aa, list):
                for i in range(len(position)):
                    # print aa, position
                    variant[position[i]] = aa[i]

                return variant
            else:
                for i in range(len(position)):
                    variant[position[i]] = aa
                return variant
        else:
            variant[position] = aa
            return variant

    def get_aa_list(self, position, bias=False, bias_type='exclude', **Kargs):
        """Helper function. Return list of AA for a given position.
        Remove from the list of potential mutations the
        current AA in order to do not produce silent mutations, if
        bias is activated, a list of AminoAcids included in bias,
         is excluded if the bias type is only, then only the aa on
         the list are considered as a option.

        Paramenters
        -----------

        postion : int
            Position to change. from: 0->x.

        bias : list
            list of aminoacids we do want to exclude or use

        bias_type : str
            how to handle the bias AA,  exclude or restrict , by default exclude

        Returns
        -------

            list
            Returns a list with the potential aminoacids




        """
        try:
            aa_alternatives_list = list(self._aa_list)
            a = self.seq[position]
            aa_alternatives_list.remove(a)
        except:
            print(sys.exc_info())
            print('{}\t{}\t{}'.format(position, self.seq, self.len_))
            sys.exit()

        if bias and (isinstance(bias, list) or isinstance(bias, tuple)):

            if bias_type == 'exclude':

                for aminoacid in bias:
                    # current postion aa check
                    if a != aminoacid:
                        aa_alternatives_list.remove(aminoacid)
            elif bias_type == 'restrict':

                aa_alternatives_list = bias
                try:
                    aa_alternatives_list.remove(a)
                except:
                    pass

            else:
                raise KeyError(
                    'bias_type unsupported: exclude or restrict are the only variables addmited')

        # print aa_alternatives_list,a,position,'aqui'
        return aa_alternatives_list

    def random_pair_mutant(self, positions, **Kargs):
        """Return a seq with random mutations in a given list of positions
        positions: list of integers 0->x

        Paramenters
        -----------

        postions : list
            Positions list of positions to change. from: 0->x.

        bias : list
            list of aminoacids we do want to use


        Returns
        -------
            str
            sequence with the modifications



        """
        changes = []
        for pos in positions:
            # print self.get_aa_list(pos)
            changes.append(random.choice(self.get_aa_list(pos, **Kargs)))
        # print changes, 'primo'
        return self.mod_template(positions, changes)

    def get_all_mutants_in(self, position, **Kargs):
        '''Return all the mutated sequences with the given position modified

        Paramenters
        -----------

        :param postion: Position to change. from: 0->x.
        :rtype: int

        :param bias: list of aminoacids we do want to use
        :rtype: list

        :yields: str -- Yields a str  with the all the mutations
        '''

        list_of_subtitutions = self.get_aa_list(position, **Kargs)

        for a in list_of_subtitutions:
            if a != self.seq[position]:
                yield self.mod_template(position, a)

    def soft_randomization(self, num_mutations=2, inrange=False, bias=False):
        """Randomly select N postions (num_mutations), and randomly
        mutanted and return seq.

        Paramenters
        -----------

        num_mutations : int
            Number of positions to mutate

        inrange : list, or tuple
            list or tuple with [start,end] coordinates inside the seq
            (by default full length)

        bias : list
            ['A','W'] aminoacids excluded

        Returns
        -------

            str Sequences with the mutations


         """
        positions = []
        if not inrange:
            start = 0
            end = self.len_ - 1
        else:
            start = inrange[0] - 1
            end = inrange[1] - 1

        peptide_positions = list(range(start, end + 1))

        for i in range(0, num_mutations):
            muta_at = (random.choice(peptide_positions))
            peptide_positions.remove(muta_at)
            positions.append(muta_at)

        random_mutation = self.random_pair_mutant(positions, bias=bias)
        random_mutation = ''.join(random_mutation)

        return random_mutation

    def frame_mutations(self, frame_size=2, inrange=False, shared_bias=False, **Kargs):
        """This method, return a sequence of corraleted mutations. Need improvement
        suport more than a pair frame, allow  specific bias for each position.

        Paramenters
        -----------

        Yields
        ------

        """

        if not inrange:
            start = 0
            end = self.len_ - 1
        else:
            start = inrange[0] - 1
            end = inrange[1] - 1

        for i in range(start, end + 1):
            if i != end:
                for aa1 in self.get_aa_list(i, bias=shared_bias, **Kargs):
                    for aa2 in self.get_aa_list(i + 1, bias=shared_bias, **Kargs):

                        yield self.mod_template([i, i + 1], [aa1, aa2])


def check_lib_integrty(file_library, master_seq, oligo_length,
                       CONSTANT_F='CAGCCTCTTCATCTGGC',
                       CONSTANT_R='GGTGGAGGATCCGGAG',
                       restriction_site=['GAATTC', 'CCCGGG', 'GTCGAC']):
    """Check if the librar complains with restrictions sites and lenght setup.

    Parameters
    ----------

    Returns
    -------



    """
    # create dict seq
    dict_seq = {}

    with open(file_library, 'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                line = line.strip()
                pep, mut, num, aaseq = line.split('_')

                seq = next(input_file)
                seq = seq.strip()
                seq = seq.upper()
    #            print seaq

                # check length
                if len(seq) != oligo_length:
                    print('ERROR Length on {} {} {} {}'.format(pep, mut, num, seq))

                # check restriction sites

                for site in restriction_site:
                    if seq.find(site) != -1:
                        print('RESTRICTION SITE FOUND {} at {} on {} {} {} {}'.format(site,
                                                                                      seq.find(
                                                                                          site),
                                                                                      pep,
                                                                                      mut,
                                                                                      num,
                                                                                      seq))

                skip = len(CONSTANT_F)
                template = translate2aa(seq[skip:len(aaseq) * 3 + skip])

                if aaseq != template:
                    print('ERROR on translation  {} {} {} {} {} {}'.format(pep, mut, num, aaseq,
                                                                           template, seq))

                if mut == 'SV':
                    if ltein.distance(master_seq, template) != 1:
                        print('ERROR  more mut than spectected in a SV {} {} {} {} {} {}'.format(pep,
                                                                                                 mut,
                                                                                                 num,
                                                                                                 master_seq,
                                                                                                 template,
                                                                                                 seq))

                if mut == 'RA':
                    if ltein.distance(master_seq, template) != 2:
                        print('ERROR  more mut than spectected in a RA {} {} {} {} {} {}'.format(pep,
                                                                                                 mut,
                                                                                                 num,
                                                                                                 master_seq,
                                                                                                 template,
                                                                                                 seq))

                if seq.find(CONSTANT_F) != 0:
                    print('ERROR with the contatn region R {} {} {} {} {} {}'.format(pep,
                                                                                     mut,
                                                                                     num,
                                                                                     aaseq,
                                                                                     template,
                                                                                     seq))

                if seq.find(CONSTANT_R) + len(CONSTANT_R) != oligo_length:
                    print('ERROR with the contatn region L {} {} {} {} {} {}'.format(pep,
                                                                                     mut,
                                                                                     num,
                                                                                     aaseq,
                                                                                     template,
                                                                                     seq))
    #                print seq.find(constan1t_L)+len(constant_L)
    #            if len(template) != len(AA):
    #                print 'ERROR with the len', pep,mut,num,aaseq,template,seq
                if template in dict_seq:
                    print('ERROR seq duplicated {} {} {} {} {} {}'.format(pep,
                                                                          mut,
                                                                          num,
                                                                          aaseq,
                                                                          template,
                                                                          seq))

                else:
                    dict_seq[template] = 1


# quick methods

# Alanine scanning Library
def Alanine_scanning(pep_library, AA='A', **Kargs):

    # aminoacids = ['I','L','V','F',' M','C','G','A','P','T','S','Y','W','Q','N',
    #             'E','D','K','R']
    # aminoacids.remove(AA)
    pep_library.generate_single_variants(bias='A', bias_type='restrict', **Kargs)

    return pep_library

# need to be improved
# Postion Scanning Library
# def position_scanning(pep_library,inrange)
# screen_library(inrage=)

# Overlapping Peptide Library
# def overlapping_pep_lib(pep_library,overlap,**Kargs):

#    pep_library.generate_inframe_variants(self,frame_size=overlap,**Kargs)

#    return pep_library

# change inframe for linked mutations


def setup_options():

    parser = argparse.ArgumentParser("""
    Generation Libraries

    Usage: %prog -t TEMPLATE_SEQ -options
    Run `%prog -h` for help information
    """)

    # actions
    # Kind of generations
    # Parameters of the generators
    # restriction ezymes
    # constant regions
    # config file?
    parser.add_option('-t', '--template', dest='template_seq',
                      default=False, help='')

    (opts, args) = parser.parse_args()

    if not opts.template_seq:
        parser.error('A template seq  not given')

    return opts, args

if __name__ == '__main__':

    # test_library_mutations()
    pass

# TODO
# implement Quick methods
# Alanine scanning
#
# Postion Scanning Library
# Random library
# scrambled Library

# Options functionalality

# import pandas as pd
# pd.read_cvs('ecoli.freq',delim_whitespace=True)
#
# LOG_FILENAME = 'default.log'
# logging.basicConfig(filename=LOG_FILENAME,
#                    level=logging.DEBUG,
#                    )
# logg = logging.getLogger(__name__)
# logg.debug('This message should go to the log file')
#
#
