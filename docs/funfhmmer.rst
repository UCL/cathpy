cathpy.funfhmmer
================

.. code:: python

    from cathpy import funfhmmer

    api = funfhmmer.Client()

    res = api.search_fasta(fasta_file='/path/to/seq.fa')

    res.funfam_scan.as_tsv()

    res.funfam_resolved_scan.as_tsv()

.. automodule:: cathpy.funfhmmer
    :members:

