.. _DockerPlaybook:

SEEK Docker Playbook
====================

To enable workflows to run seamlessly across workstations without
lengthy and buggy installations, we utilize Docker containers that house
specific software that are utilized by different rules in ``snakemake`` workflows.

See `docker hub <https://hub.docker.com/orgs/neuroseek/repositories>`_ for
full list of Docker containers we maintain.

Updating Docker Containers
--------------------------
*Before starting new code*, we highly recommend opening an issue on `GitHub <https://github.com/ncsl/seek>`_ to discuss potential changes.
Please, also see the `Contributing Guide <https://github.com/ncsl/seek/CONTRIBUTING.md>`_.

To update Docker containers, one needs to check:

1. Update the ``Dockerfile`` recipe inside ``/dockerfiles/`` of the repository. For example,
``Dockerfile.acpcdetect`` hosts the Docker recipe for running ``acpcdetect``.
2. Are any of the workflows affected? To check this, run unit tests.
3. Do the deployment recipes need to be modified? Check the `Makefile <https://github.com/ncsl/seek/Makefile>`_
to see if any recipes or versions of containers need to be updated.

Building and Pushing Docker Containers
--------------------------------------

We include a ``Makefile`` recipe for building and pushing each Docker container. For example

    .. code-block:: bash

        $ make build-acpcdetect

will build the ``acpcdetect`` Docker container. Then to push it to ``neuroseek`` Docker Hub,
one must have push privileges (for developers)

    .. code-block:: bash

        $ make push-acpcdetect

To build all images at once that are used by seek, we use a ``docker-compose.yml`` file
to build all the images

    .. code-block:: bash

        $ make build-all

        # build the composition in `docker-compose.yml`
        $ docker-compose up --build

        # run the container
        $ docker-compose up

To pull the images all at once instead of building them locally

    .. code-block:: bash

        $ make pull-all

To push all the images that you have tested and built locally

    .. code-block:: bash

        $ make push-all

Testing Docker Containers
-------------------------
If any Docker containers need to be updated, then the workflows should be tested. See the
`testing guide <testing>_`.
