import os
import pytest


from ngskit.demultiplex_reads import *

@pytest.mark.parametrize("seq, target,cutoff,expected", [
    ('ATGATGAAAA', 'ATGATGAAAA', 0, True),
    (b'ATGATGAAAA', 'ATGATGAAAG', 0, False), 
    (b'ATGATGAAAA', 'ATGATGAAAG', 1, True), 
    (b'ATGATGAAAA', 'ATGATGAGGG', 3, True),
    ('ATGATGAAAA', 'ATGATGAGGG', 3, True),
    ('ATGATGAAAA', 'ATGATGAGGG', 2, False),
])
def test_eval(seq, target,expected):
    assert eval(match(seq, target)) == expected