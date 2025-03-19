Contributing
==============

Thank you for contributing!


Reporting Issues
------------------

To report a bug or request a feature please open an 
[issue](https://github.com/seshadrilab/tcrconvertr/issues).


Contributing Code
-------------------

**1. Install suggested dependencies**

Install the testing and documentation dependencies listed in ``pyproject.toml``:

.. code-block:: console

   $ pip install .[dev]
   $ pip install .[docs]

**2. Fork the repo and make changes**

- Fork the repository and create a branch off of `main`.
- Ensure changes are covered by tests.
- Update the documentation as needed.

**3. Run checks**

GitHub Actions will perform linting and run package checks and tests when you 
push changes. You can also check your code ahead of time:

.. code-block:: console

   # Tests
   $ pytest

   # Run code examples
   $ python -m doctest <changed_script.py>

   # Linting, change format of files to match style
   $ ruff format

**4. When ready, open a pull request (PR)**

- Include a clear description of the changes.
- Reference any related issues.
- Make sure all checks pass.
