import os
import pytest

from NGSKit import oligolib_generator as olg


def test_lib_default_creation():
    """Create default library"""
    test_lib = olg.peptide_library('AAAA')

    assert test_lib.lib_name == 'AAAA'
    assert test_lib.specie == 'human'
    assert test_lib.include_template_in
    assert len(test_lib) == 1
    assert test_lib._aalibrary['AAAA'] == 'AAAA_OO_0'


def test_lib_custom_creation():

    # create a new library abd check if user set values works
    test_lib = olg.peptide_library('AAA', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    assert test_lib.lib_name != 'AAA'
    assert test_lib.specie != 'human'
    assert test_lib.include_template_in == False
    assert test_lib.CONSTANT_R == 'AGT'
    assert test_lib.CONSTANT_L == 'XXX'
    assert len(test_lib) == 0


def test_lib_creation_from_a_file():

    # test_lib to load designs and add restriction enzymes
    test_lib = olg.peptide_library('AAA', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')
    test_lib.load_designs('./test_file.inp')
    assert len(test_lib) == 2
    assert test_lib._aalibrary['AAC']
    assert test_lib._aalibrary['AAI']


def test_add_restriction_enzymes():

    test_lib = olg.peptide_library('AAA', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')
    test_lib.load_designs('./test_file.inp')

    test_lib.add_restriction_enzyme(['GCTGCC', 'GCAGCT'])
    test_lib.add_restriction_enzyme({'free': 'GCAGCG'})
    assert len(test_lib.restriction_enzyme) == 6


def test_write_library():

    # create a lib from a file
    test_lib = olg.peptide_library('AAA', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')
    test_lib.load_designs('./test_file.inp')
    # add restriction enzymes
    test_lib.add_restriction_enzyme(['GCTGCC', 'GCAGCT'])
    test_lib.add_restriction_enzyme({'free': 'GCAGCG'})
    # test_lib write functions
    test_lib.write(file_name='wtest0.fasta')
    # test different restriction enzymes
    test_lib.add_restriction_enzyme(['GCCGCG', 'GCTGCT'])
    test_lib.write(file_name='wtest1.fasta')

    test_lib.add_restriction_enzyme(['GCTGCG'])
    test_lib.write(file_name='wtest2.fasta')

    test_lib.add_restriction_enzyme(['GCTGCA'])
    test_lib.write(file_name='wtest3.fasta')

    test_lib.add_restriction_enzyme(['GCGGCG'])
    test_lib.write(file_name='wtest4.fasta')

    # assert files
    for i in range(5):
        assert os.path.isfile('wtest{}.fasta'.format(i))

    # rmove test files
    for i in range(5):
        os.remove('wtest{}.fasta'.format(i))


def test_change_template():

    test_lib = olg.peptide_library('AAAA')
    test_lib.load_designs('./test_file.inp')
    test_lib.change_template('CCC')
    assert test_lib.template.seq == 'CCC'


def test_simple_var_generation():
    # test_lib generation of all simple variants
    # by default option
    test_lib = olg.peptide_library('CCC', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_single_variants()
    test_lib.write(file_name='wtestvar.fasta')
    assert len(test_lib) == 57
    assert os.path.isfile('wtestvar.fasta')
    os.remove('wtestvar.fasta')


def test_bias_var_generation():
    # test_lib generation of all variant with exclusion variant
    test_lib = olg.peptide_library('CCC', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_single_variants(bias=['A'])
    assert len(test_lib) == 54


def test_restrict_double_bias_var_generation():
    # test_lib generation of all variant with restriction variants
    test_lib = olg.peptide_library('CCC', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_single_variants(bias=['G', 'C'], bias_type='restrict')
    assert len(test_lib) == 3


def test_double_bias_var_generation():
    # test_lib generation of all variant with exclusion variant more than one
    test_lib = olg.peptide_library('CCC', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_single_variants(bias=['A', 'G'])
    assert len(test_lib) == 51

    # check if the library is free of Ala and Gly
    for seq in test_lib._aalibrary.keys():
        assert not 'A' in seq
        # raise 'error in bias'
        assert not 'G' in seq
        # raise 'error in bias'


def test_writing_options():
    # Test writing
    test_lib = olg.peptide_library('CCC',  include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_single_variants(bias=['A', 'G'])
    # extend sequence to 25
    test_lib.write(file_name='test-lib.fasta', extend=25)

    olg.check_lib_integrty('test-lib.fasta', 'CCC', 25, CONSTANT_R='AGT',
                           CONSTANT_L='XXX', restriction_site=test_lib.restriction_enzyme.values())

    os.remove('test-lib.fasta')


def test_random_trim():
    # Test random removing
    test_lib = olg.peptide_library('CCC', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_single_variants(bias=['A', 'G'])
    test_lib.random_remove(21)
    assert len(test_lib) == 30


def test_radom_generation():
    # test_lib random generation
    test_lib = olg.peptide_library('CCC',  include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')

    test_lib.generate_random_variants(how_many=10, mutant_kind=[2])
    test_lib.write(add_stop_end=False)
    olg.check_lib_integrty('test-lib.fasta', 'CCC', 15, CONSTANT_R='AGT',
                           CONSTANT_L='XXX', restriction_site=test_lib.restriction_enzyme.values())

    # test_lib random generation in range
    test_lib._aalibrary = dict()
    test_lib.generate_random_variants(how_many=10, mutant_kind=[2], inrange=[2, 3])

    # test_lib write without extend and stop codons
    # test_lib inframe
    # All sequence are the same at the beging
    test_lib.write(add_stop_end=False)
    assert os.path.isfile('test-lib.fasta')
    os.remove('test-lib.fasta')

    test_lib._aalibrary = dict()

    test_lib.generate_inframe_variants(frame_size=2)
    test_lib.write(add_stop_end=False)

    olg.check_lib_integrty('test-lib.fasta', 'CCC', 15, CONSTANT_R='AGT',
                           CONSTANT_L='XXX', restriction_site=test_lib.restriction_enzyme.values())

    os.remove('test-lib.fasta')



def test_load_json():
    test_lib = olg.peptide_library('CCC', include_template=False, CONSTANT_R='AGT',
                                   CONSTANT_L='XXX', lib_name='test-lib', specie='E.coli')


    # test_lib json generation permutations

    test_lib.permutations_json('test_permutations_cfg.json')
    assert len(test_lib._aalibrary) == 6
    test_lib.write()
    assert os.path.isfile('test-lib.fasta')
    os.remove('test-lib.fasta')
