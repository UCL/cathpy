cathpy.core.align
=================

.. code:: python

    from cathpy.core.align import Align

    aln = Align.from_stockholm('/path/to/align.sto')

    aln.count_sequences
    # 75

    seq = aln.find_seq_by_accession('1cukA01')
    seq = aln.find_seq_by_id('1cukA01/1-151')


.. automodule:: cathpy.core.align
    :members:


