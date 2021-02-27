PYTHON ?= python
PYTESTS ?= pytest
CODESPELL_SKIPS ?= "seek/pipeline/*,*_old/*,doc/auto_*,*.fif,*.eve,*.gz,*.tgz,*.zip,*.mat,*.stc,*.label,*.w,*.bz2,*.annot,*.sulc,*.log,*.local-copy,*.orig_avg,*.inflated_avg,*.gii,*.pyc,*.doctree,*.pickle,*.inv,*.png,*.edf,*.touch,*.thickness,*.nofix,*.volume,*.defect_borders,*.mgh,lh.*,rh.*,COR-*,FreeSurferColorLUT.txt,*.examples,.xdebug_mris_calc,bad.segments,BadChannels,*.hist,empty_file,*.orig,*.js,*.map,*.ipynb,searchindex.dat,install_mne_c.rst,plot_*.rst,*.rst.txt,c_EULA.rst*,*.html,gdf_encodes.txt,*.svg"
CODESPELL_DIRS ?= seek/ workflow/ doc/ tutorials/

all: clean inplace test

# variables
name := seek
version := 1.0.1
dockerhub := neuroseek

# docker containers
acpcdetect_version := 2.0
freesurfer7-with-mrtrix3_version := 1.2

############################## MAIN COMMANDS #########################

snakemake-all: recon coregistration prep_viz

# snakemake command line optional arguments
cores := 1
sing-args := "--bind ~/hdd/epilepsy_bids/,/home/adam2392/Documents/seek/"

recon:
	cd workflow/recon && \
	snakemake --cores $(cores) --use-singularity --singularity-prefix ../.singularity/ --singularity-args $(sing-args);

prep-localization:
	cd workflow/prep_localization && \
	snakemake --cores $(cores) --use-singularity --singularity-prefix ../.singularity/ --singularity-args $(sing-args);

coregistration:
	cd workflow/coregistration && \
	snakemake --cores $(cores) --use-singularity --singularity-prefix ../.singularity/ --singularity-args $(sing-args);

prep-viz:
#	@read -p "Enter full absolute path to 'seek' repository:" path;
#module_dir=./modules/$$module;
	cd workflow/prep_vizengine && \
	snakemake --cores $(cores) --use-singularity --singularity-prefix ../.singularity/ --singularity-args $(sing-args);

############################## DOCKER #########################
build-acpc:
	cp ./.data/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz ./dockerfiles/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz
	docker build --rm -f ./dockerfiles/Dockerfile.acpcdetect -t $(dockerhub)/acpcdetect:$(acpcdetect_version)  ./dockerfiles
	rm ./dockerfiles/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz

build-freesurfer:
	docker build --rm -f ./dockerfiles/Dockerfile.freesurfer-with-mrtrix3 -t $(dockerhub)/freesurfer7-with-mrtrix3:$(freesurfer7-with-mrtrix3_version)  ./dockerfiles

build-seek:
	cp Pipfile ./dockerfiles/Pipfile
	docker build --rm -f ./dockerfiles/Dockerfile.seek -t $(dockerhub)/seek:$(version)  ./dockerfiles
	rm ./dockerfiles/Pipfile

push-acpc:
	docker push $(dockerhub)/acpcdetect:$(acpcdetect_version)

push-freesurfer:
	docker push $(dockerhub)/freesurfer7-with-mrtrix3:$(freesurfer7-with-mrtrix3_version)

push-seek:
	docker push $(dockerhub)/seek:$(version)

pull-all:
	docker pull $(dockerhub)/acpcdetect:$(acpcdetect_version)
	docker pull $(dockerhub)/blender:$(blender_version)
	docker pull $(dockerhub)/freesurfer7-with-mrtrix3:$(freesurfer7-with-mrtrix3_version)
	docker pull docker://cbinyu/fsl6-core

############################## UTILITY FOR SNAKEMAKE #########################
outputpath := "./doc/_static"

init:
    export SEEKHOME=$(shell pwd)

create_dags:
	snakemake --snakefile ./workflow/recon_workflow/Snakefile --forceall --dag | dot -Tpdf > $(outputpath)/recon_workflow.pdf;
	snakemake --snakefile ./workflow/prep_localization_workflow/Snakefile --forceall --dag | dot -Tpdf > $(outputpath)/prep_localization_workflow.pdf;
	snakemake --snakefile ./workflow/coregistration_workflow/Snakefile --forceall --dag | dot -Tpdf > $(outputpath)/coregistration_workflow.pdf;
	snakemake --snakefile ./workflow/prep_vizengine_workflow/Snakefile --forceall --dag | dot -Tpdf > $(outputpath)/prep_viz_workflow.pdf;

############################## UTILITY FOR PYTHON #########################
clean-build:
	rm -rf _build

clean-ctags:
	rm -f tags

clean-cache:
	find . -name "__pychache__" | xargs rm -rf

clean: clean-build clean-ctags clean-cache

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

snakelint:
	cd workflow/recon/ && snakemake --lint;
	cd workflow/prep_localization/ && snakemake --lint;
	cd workflow/coregistration/ && snakemake --lint;
	cd workflow/prep_vizengine/ && snakemake --lint;

check:
	@$(MAKE) -k black pydocstyle codespell-error

