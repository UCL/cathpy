cathpy.version
==============

.. code:: python

    from cathpy.version import CathVersion

    cv = CathVersion("v4.2") # or "v4_2_0", "current"

    cv.dirname
    # "4_2_0"

    cv.pg_dbname
    # "cathdb_v4_2_0"

    cv.is_current
    # False

.. automodule:: cathpy.version
    :members:

