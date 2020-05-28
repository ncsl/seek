:Overview: Preps the Reconstruction and BIDS Layout
:Input: T1 MRI Dicoms, CT Dicoms (optional)
:Output:
    - count.txt

Usage
~~~~~~~

::

    sequana init pipeline_count
    snakemake -s pipeline_count.rules -f

Requirements
~~~~~~~~~~~~~~

Here you should list the dependencies, which should match the file
requirements.txt in ./sequana_pipelines/count/

.. image:: https://raw.githubusercontent.com/sequana/sequana_count/master/sequana_pipelines/count/dag.png

Details
~~~~~~~~~~~~~

Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

count rule
^^^^^^^^^^^
.. snakemakerule:: pipeline_count