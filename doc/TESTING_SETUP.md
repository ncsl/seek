# Guide to setup Testing Environment

# Testing Scipts

    black neuroimg/*
    black tests/*
    black --check neuroimg/
    pytest tests/ 
    coverage-badge -f -o coverage.svg

# Updating Packages

    conda update --all
    conda env export > environment.yaml

# Documentation 
# TODO: setup guides for running autodoc

    sphinx-quickstart