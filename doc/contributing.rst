.. _Contribute:

Contributing to SEEK
====================

(adopted from scikit-learn)

The latest contributing guide is available in the repository at
`doc/contributing.rst`, or online at:

There are many ways to contribute to NeuroimgPipe, with the most common ones
being contribution of code or documentation to the project. Improving the
documentation is no less important than improving the pipeline itself. If you
find a typo in the documentation, or have made improvements, do not hesitate to
submit a GitHub pull request. Documentation can be found under the
[doc/](https://github.com/ncsl/seek/tree/master/doc) directory.

But there are many other ways to help. In particular answering queries on the
[issue tracker](https://github.com/ncsl/seek/issues), and
investigating bugs are very valuable contributions that decrease the burden on 
the project maintainers.

Another way to contribute is to report issues you're facing, and give a "thumbs
up" on issues that others reported and that are relevant to you. It also helps
us if you spread the word: reference the project from your blog and articles,
link to it from your website, or simply star it in GitHub to say "I use it".

Another way to contribute is specifically to make additional pipelines that improve 
the accuracy of contact localization for iEEG data using the T1 and CT images.

Code of Conduct
---------------

We abide by the principles of openness, respect, and consideration of others
of the Python Software Foundation: https://www.python.org/psf/codeofconduct/.

Code Guidelines
----------------

*Before starting new code*, we highly recommend opening an issue on `GitHub <https://github.com/ncsl/seek>`_ to discuss potential changes.

* Please use standard `black <https://black.readthedocs.io/en/stable/>`_ Python style guidelines. To test that your code complies with those, you can run:

  .. code-block:: bash

     $ black --check seek/
     $ make check

  In addition to `black`, `make check` command runs `pydocstyle` and `codespell`

* Use `NumPy style <https://numpydoc.readthedocs.io/en/latest/format.html>`_ for docstrings. Follow existing examples for simplest guidance.

* New functionality must be **validated** against sample datasets.

* Changes must be accompanied by **updated documentation** and examples.

* After making changes, **ensure all tests pass**. This can be done by running:

  .. code-block:: bash

     $ pytest --doctest-modules

Snakemake Rule Guidelines
-------------------------

*Before starting a new rule*, we highly recommend opening an issue on `GitHub <https://github.com/ncsl/seek>`_ to discuss potential changes.

* For `snakemake` rules, we use https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html

Snakemake utilizes a set of configuration files for the user to define their dataset. This is done by editing the `subjects.tsv` file
in the `pipeline/config/` directory. Configuration files are checked against the designed schemas in `pipeline/schemas` directory.
Each corresponding sub-workflow is documented in a `.smk` Snakemake file, which can define things like `Prep`, which prepares
raw dicom files from scratch to run through the SEEK workflow.
