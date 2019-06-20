cath-align-scorecons
====================

.. code:: python

  usage: cath-align-scorecons [-h] --in IN_FILE [--out OUT_FILE] --format
                              {sto,fasta} [--replace] [--verbose]

  Update the scorecons data for a STOCKHOLM alignment

  optional arguments:
    -h, --help            show this help message and exit
    --in IN_FILE, -i IN_FILE
                          input alignment file (default: None)
    --out OUT_FILE, -o OUT_FILE
                          output STOCKHOLM alignment file (default: None)
    --format {sto,fasta}, -f {sto,fasta}
                          alignment format of input file (default: sto)
    --replace             overwrite the input alignment file with updated
                          scorecons (default: False)
    --verbose, -v         more verbose logging (default: 0)
