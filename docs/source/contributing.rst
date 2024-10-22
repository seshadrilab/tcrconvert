Contributing
==============

To contribute, please submit a pull request on our GitHub `page <https://github.com/seshadrilab/tcrconvert/pulls>`_.


Testing
---------

It's essential that you write tests for any new code that performs actions (not documentation changes). 
You should also run the current test suite to ensure your code doesn't break anything.

To get started, first install the necessary testing and documentation dependencies listed in ``requirements-dev.txt``:

.. code-block:: console

   $ pip install -r requirements-dev.txt

Then you can run tests using ``pytest`` from the root of the GitHub repository:

.. code-block:: console

   $ pytest

Additionally, the ``tests`` and ``codecov`` GitHub Actions automatically run the current test suite and 
check for test coverage on every push and pull request. 


Linting
---------

We recommend `Ruff <https://docs.astral.sh/ruff/>`_ for linting. You can install it from PyPi:

.. code-block:: console

   $ pip install ruff

The ``ruff`` GitHub Action will automatically perform linting on every push and pull request.