cathpy.core.models
==================

Provides access to classes that representing general entities such as amino acids, db identifiers, etc.

.. code:: python

    from cathpy.core.models import (
        AminoAcid,
        AminoAcids, 
        CathID,
        ClusterFile,
        FunfamID, )

    aa = AminoAcids.get_by_id('A')

    aa.one                      # 'A'
    aa.three                    # 'ala'
    aa.word                     # 'alanine'

    AminoAcids.is_valid_aa('Z') # False

    cathid = CathID("1.10.8.10.1")

    cathid.sfam_id              # '1.10.8.10'
    cathid.depth                # 5
    cathid.cath_id_to_depth(3)  # '1.10.8'
    
    funfam_file = ClusterFile("/path/to/1.10.8.10-ff-1234.reduced.sto")

    funfam_file.path            # '/path/to/'
    funfam_file.sfam_id         # '1.10.8.10'
    funfam_file.cluster_type    # 'ff'
    funfam_file.cluster_num     # 1234
    funfam_file.desc            # '.reduced'
    funfam_file.suffix          # '.sto'

.. automodule:: cathpy.core.models
    :members:

