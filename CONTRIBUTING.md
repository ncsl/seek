# Contributions

Contributions are welcome in the form of pull requests.

Once the implementation of a piece of functionality is considered to be bug
free and properly documented (both API docs and an example script),
it can be incorporated into the master branch.

To help developing `seek`, you will need a few adjustments to your
installation as shown below.

## Running tests

### (Optional) Install Docker
To run workflows, it is recommended to use Docker. Install Docker at: https://docs.docker.com/get-docker/.

### Install development version of seek
First, you should [fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the `seek` repository. 
Then, clone the fork and install it in.

TODO insert how to install and setup, or link to instructions.

### Install Python packages required to run tests
Install the following packages for testing purposes, plus all optonal MNE-BIDS
dependencies to ensure you will be able to run all tests.

    $ pip install flake8 pytest pytest-cov

Running code linter:

    

## Building the documentation

TODO: @Christopher Coogan, can you help fill this in?

The documentation can be built using sphinx. For that, please additionally
install the following:

    $ pip install matplotlib nilearn sphinx numpydoc sphinx-gallery sphinx_bootstrap_theme pillow

To build the documentation locally, one can run:

    $ cd doc/
    $ make html

or

    $ make html-noplot
    
if you don't want to run the examples to build the documentation. This will result in a faster build but produce no plots in the examples.

### Creating a DAG for each workflow
``Snakemake`` has the nice property of implicitly building a directed acyclic graph (DAG) of 
your anatomical analysis. To visualize this DAG for your subjects, one can make use of the 
snakemake command:

    $ snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > ./recon_workflow.pdf;

We package all these DAG creation commands into a ``Makefile`` recipe, which can be ran as:

    $ make create_dags

## BIDS-Validation
To robustly apply seek workflows and reconstruction visualiztion, we rely on the BIDS specification 
for storing data. One can use the `bids-validator` to verify that a dataset is BIDS-compliant.

### Install the BIDS validator
Finally, it is necessary to install the
[BIDS validator](https://github.com/bids-standard/bids-validator). The outputs
of MNE-BIDS are run through the BIDS validator to check if the conversion
worked properly and produced datasets that conforms to the BIDS specifications.

You will need the `command line version` of the validator.

#### Global (system-wide) installation
- First, install [Node.js](https://nodejs.org/en/).
- For installing the **stable** version of `bids-validator`, please follow the
instructions as detailed in the README of the bids-validator repository.
- For installing the **development** version of `bids-validator`, see [here](https://github.com/bids-standard/bids-validator/blob/master/CONTRIBUTING.md#using-the-development-version-of-bids-validator).

Test your installation by running:

    $ bids-validator --version

#### Local (per-user) development installation

Install [Node.js](https://nodejs.org/en/). If you're use `conda`, you can
install the `nodejs` package from `conda-forge` by running
`conda install -c conda-forge nodejs`.

Then, retrieve the validator and install all its dependencies via `npm`.

    $ git clone git@github.com:bids-standard bids-validator.git
    $ cd bids-validator/bids-validator
    $ npm i

Test your installation by running:

    $ ./bin/bids-validator --version
