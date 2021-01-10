:orphan:

.. _whats_new:


What's new?
===========

Here we list a changelog of SEEK.

.. contents:: Contents
   :local:
   :depth: 3

.. currentmodule:: seek

.. _current:

Current
-------

.. _changes_0_1:

Version 0.1
-----------

Changelog
~~~~~~~~~

- Refactored semi-automated algorithm for localizing contacts on CT img, in :code:`seek/localizae_contacts/electrode_clustering` by `Chester Huynh`_ (:gh:`16`)
- Pipeline for group comparisons with the MNI152 atlas at :code:`seek/pipeline/07-group_analysis` by `Adam Li`_ (:gh:`25`)
- Add full set of installation, and usage instructions for seek pipeline, and example DAG pdfs to documentation by `Adam Li`_ (:gh:`26`)

Bug
~~~


API
~~~
- Added a Makefile recipe for running all the related SEEK workflows, by `Adam Li`_ (:gh:`34`)


Authors
~~~~~~~

People who contributed to this release (in alphabetical order):

* Adam Li
* Chester Huynh
* Christopher Coogan

:doc:`Find out what was new in previous releases <whats_new_previous_releases>`

.. include:: authors.rst