PYTHON ?= python
PYTESTS ?= pytest
CODESPELL_SKIPS ?= "seek/pipeline/*,*_old/*,doc/auto_*,*.fif,*.eve,*.gz,*.tgz,*.zip,*.mat,*.stc,*.label,*.w,*.bz2,*.annot,*.sulc,*.log,*.local-copy,*.orig_avg,*.inflated_avg,*.gii,*.pyc,*.doctree,*.pickle,*.inv,*.png,*.edf,*.touch,*.thickness,*.nofix,*.volume,*.defect_borders,*.mgh,lh.*,rh.*,COR-*,FreeSurferColorLUT.txt,*.examples,.xdebug_mris_calc,bad.segments,BadChannels,*.hist,empty_file,*.orig,*.js,*.map,*.ipynb,searchindex.dat,install_mne_c.rst,plot_*.rst,*.rst.txt,c_EULA.rst*,*.html,gdf_encodes.txt,*.svg"
CODESPELL_DIRS ?= seek/ workflow/ doc/ tutorials/

all: clean inplace test

# variables
name := seek
version := 1.0.0
dockerhub := neuroseek

# docker containers
blender_version := 2.82
acpcdetect_version := 2.0
freesurfer7-with-mrtrix3_version := 1.2

############################## MAIN COMMANDS #########################
snakemake_all: recon coregistration prep_viz

recon:
	snakemake --cores 1 --use-singularity --singularity-args "--bind ~/hdd/epilepsy_bids/"

coregistration:
	snakemake --cores 1 --use-singularity --singularity-args "--bind ~/hdd/epilepsy_bids/"

prep_viz:
	snakemake --cores 1 --use-singularity --singularity-args "--bind ~/hdd/epilepsy_bids/"


############################## DOCKER #########################
build:
	@docker-compose build;

build-acpc:
	docker build --rm -f ./dockerfiles/Dockerfile.acpcdetect -t $(dockerhub)/acpcdetect:$(acpcdetect_version)  ./dockerfiles

build-blender:
	docker build --rm -f ./dockerfiles/Dockerfile.meshgenerator -t $(dockerhub)/blender:$(blender_version)  ./dockerfiles

build-freesurfer:
	docker build --rm -f ./dockerfiles/Dockerfile.freesurfer-with-mrtrix3 -t $(dockerhub)/freesurfer7-with-mrtrix3:$(freesurfer7-with-mrtrix3_version)  ./dockerfiles

push-acpc:
	docker push $(dockerhub)/acpcdetect:$(acpcdetect_version)

push-blender:
	docker push $(dockerhub)/blender:$(blender_version)

push-freesurfer:
	docker push $(dockerhub)/freesurfer7-with-mrtrix3:$(freesurfer7-with-mrtrix3_version)

pull-all:
	docker pull $(dockerhub)/acpcdetect:$(acpcdetect_version)
	docker pull $(dockerhub)/blender:$(blender_version)
	docker pull $(dockerhub)/freesurfer7-with-mrtrix3:$(freesurfer7-with-mrtrix3_version)
	docker pull docker://cbinyu/fsl6-core

############################## UTILITY FOR SNAKEMAKE #########################
init:
	pipenv shell
    export SEEKHOME = $(shell pwd)

create_dags:
	snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > ./doc/_static/recon_workflow.pdf;
	snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > ./doc/_static/coregistration_workflow.pdf;
	snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > ./doc/_static/prep_viz_workflow.pdf;

############################## UTILITY FOR PYTHON #########################
clean-pyc:
	find . -name "*.pyc" | xargs rm -f
	find . -name "*.DS_Store" | xargs rm -f

clean-so:
	find . -name "*.so" | xargs rm -f
	find . -name "*.pyd" | xargs rm -f

clean-build:
	rm -rf _build

clean-ctags:
	rm -f tags

clean-cache:
	find . -name "__pychache__" | xargs rm -rf

clean: clean-build clean-pyc clean-so clean-ctags clean-cache

codespell:  # running manually
	@codespell -w -i 3 -q 3 -S $(CODESPELL_SKIPS) --ignore-words=ignore_words.txt $(CODESPELL_DIRS)

codespell-error:  # running on travis
	@echo "Running code-spell check"
	@codespell -i 0 -q 7 -S $(CODESPELL_SKIPS) --ignore-words=ignore_words.txt $(CODESPELL_DIRS)

inplace:
	$(PYTHON) setup.py install

install-tests:
	${PYTHON} -m pip install pytest black check-manifest pytest-cov pydocstyle sphinx sphinx-gallery

test: inplace check-manifest
	rm -f .coverage
	$(PYTESTS) ./
	cd seek/pipeline/01-prep/

test-doc:
	$(PYTESTS) --doctest-modules --doctest-ignore-import-errors

build-doc:
	cd doc; make clean
	cd doc; make html

pydocstyle:
	@echo "Running pydocstyle"
	@pydocstyle

pycodestyle:
	@echo "Running pycodestyle"
	@pycodestyle

check-manifest:
	check-manifest --ignore .circleci*,docs,.DS_Store,annonymize

black:
	@if command -v black > /dev/null; then \
		echo "Running black"; \
		black --check seek/; \
		black seek/; \
		black --check workflow/; \
		black workflow/; \
	else \
		echo "black not found, please install it!"; \
		exit 1; \
	fi;
	@echo "black passed"

check:
	@$(MAKE) -k black pydocstyle codespell-error

