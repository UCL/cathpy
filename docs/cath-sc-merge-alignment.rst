cath-sc-merge-alignment
=======================

.. code:: python

  usage: cath-sc-merge-alignment [-h] --cv CATH_VERSION --in SC_FILE --out
                                OUT_STO [--ff_dir FF_DIR] [--ff_tmpl FF_TMPL]
                                [--verbose]

  Merge FunFams using a structure-based alignment

  optional arguments:
    -h, --help         show this help message and exit
    --cv CATH_VERSION  cath version
    --in SC_FILE       input reference structure-based alignment (FASTA)
    --out OUT_STO      output merged alignment file (STOCKHOLM)
    --ff_dir FF_DIR    directory in which funfam alignments can be found
                      (default: <sc_dir>)
    --ff_tmpl FF_TMPL  template used to help locate funfam file (default:
                      __SFAM__-ff-__FF_NUM__.reduced.sto)
    --verbose, -v      more verbose logging
    