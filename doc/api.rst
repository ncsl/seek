:orphan:

.. _api_documentation:

=================
API Documentation
=================

Here we list the Application Programming Interface (API) for SEEK. This is the reference for classes (``CamelCase`` names) and functions
(``underscore_case`` names) of SEEK, grouped thematically by analysis stage.

.. contents:: Contents
   :local:
   :depth: 1

seek Core API (:py:mod:`seek`)
=====================================

.. currentmodule:: seek

.. autosummary::
   :toctree: generated/

    write_eztrack_bids
    label_electrode_contacts
    identify_electrode_clusters
    bids_validate
    convert_img_to_bids
    PatientBidsRoot