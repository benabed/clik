Interfacing the library with python
===================================

The module clik contains the wrapper to the clik c library.
It contains only one object called ``clik``, which is initialized with a string containing the path to a likelihood file.

.. code-block:: python

    import clik
    
    clikid = clik.clik("clikidfile")
    

The ``has_cl``, ``lmax`` and parameter names array (see :ref:`querying`) can be queried by simpliy reading the ``has_cl``, ``lmax`` and ``extra_parameter_names`` attributes of the object

.. code-block:: python

    has_cl = clikid.has_cl
    print has_cl
    
A log likelihood is computed by calling the object with a list-like object (``tuple``, ``list`` of ``numpy.ndarray`` objects) containing the vector of parameters as described in :ref:`querying`.

.. code-block:: python

    loglkl = clikid(cl_and_pars)

The file ``click_example_py.py`` gives a simple example of the use of the python API. It is compiled and installed as :program:`clik_example_py`.
