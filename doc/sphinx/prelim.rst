Preliminaries
=============

File Formats: XML and HDF5
--------------------------

XML
````

HDF5
````

Science: Units, Statistical Mechanics
-------------------------------------

Testing: Unit Tests, TDD, System Integration Tests
--------------------------------------------------

Unit Tests
``````````

Unit testing uses GoogleTest_. The unit tests are in the
`pi-qmc/unit-test/` subdirectory, which mirrors the
structure of the `pi-qmc/src` directory.

Each unit test should execute in a few miliseconds, so that
the entire suite can be run in a few seconds.

Right now the unit tests are only included in the cmake build.

.. _GoogleTest: http://code.google.com/p/googletest/

System Integration Tests
````````````````````````

System integration tests are run with python scripts. 
We use the python unittest_ module to organize the test cases.
These system tests can be run using nosetests_, like

    nosetests -v --rednose
   
System integration tests run the `pi-qmc` executable on 
real test systems, and can take a few minutes to run.

.. _unittest: http://docs.python.org/2/library/unittest.html
.. _nosetests: https://nose.readthedocs.org/en/latest/

Parallel Computing
------------------

MPI
```

OpenMP
``````
