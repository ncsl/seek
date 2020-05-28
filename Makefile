

PYTHON ?= python
PYTESTS ?= pytest
CODESPELL_SKIPS ?= "doc/auto_*,*.fif,*.eve,*.gz,*.tgz,*.zip,*.mat,*.stc,*.label,*.w,*.bz2,*.annot,*.sulc,*.log,*.local-copy,*.orig_avg,*.inflated_avg,*.gii,*.pyc,*.doctree,*.pickle,*.inv,*.png,*.edf,*.touch,*.thickness,*.nofix,*.volume,*.defect_borders,*.mgh,lh.*,rh.*,COR-*,FreeSurferColorLUT.txt,*.examples,.xdebug_mris_calc,bad.segments,BadChannels,*.hist,empty_file,*.orig,*.js,*.map,*.ipynb,searchindex.dat,install_mne_c.rst,plot_*.rst,*.rst.txt,c_EULA.rst*,*.html,gdf_encodes.txt,*.svg"
CODESPELL_DIRS ?= eztrack/ doc/ tutorials/ examples/ tests/

LOCAL_RESULTSDIR=/home/adam2392/hdd/derivatives/
EXTERNAL_RESULTSDIR=/media/adam2392/SEAGATE_HDD/derivatives/
LOCALDIR=/home/adam2392/hdd2/data/
EXTERNALDIR=/media/adam2392/SEAGATE_HDD/data/

MARCCDIR=ali39@jhu.edu@gateway2.marcc.jhu.edu:/home-1/ali39@jhu.edu/data/epilepsy_bids/
MARCC_USER=ali39@jhu.edu
ssh 							:= ssh $(port)
remote	          				:= $(MARCC_USER)@gateway2.marcc.jhu.edu

all: clean inplace test

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
    ${PYTHON} -m pip install --upgrade --no-deps https://api.github.com/repos/mne-tools/mne-python/zipball/master
    ${PYTHON} -m pip install --upgrade https://api.github.com/repos/mne-tools/mne-bids/zipball/master

test: inplace check-manifest
	rm -f .coverage
	$(PYTESTS) ./
	cd seek/pipeline/01-prep/

test-doc:
	$(PYTESTS) --doctest-modules --doctest-ignore-import-errors

test-coverage:
	rm -rf coverage .coverage
	$(PYTESTS) --cov=./ --cov-report html:coverage

trailing-spaces:
	find . -name "*.py" | xargs perl -pi -e 's/[ \t]*$$//'

build-doc:
	cd doc; make clean
	cd doc; make html

pydocstyle:
	@echo "Running pydocstyle"
	@pydocstyle

pycodestyle:
	@echo "Running pycodestyle"
	@pycodestyle

push-marcc:
	rsync -aP $(LOCALDIR) $(MARCCDIR);

push-external:
	rsync -aP $(LOCALDIR) $(EXTERNALDIR) --exclude='*/tempdir/';
	rsync -aP $(LOCAL_RESULTSDIR) $(EXTERNAL_RESULTSDIR);

pull-marcc:
	rsync -aP $(MARCCDIR) $(LOCALDIR);

pull-external:
	rsync -aP $(EXTERNALDIR) $(LOCALDIR) --exclude='*/tempdir/';

check-manifest:
	check-manifest --ignore .circleci*,docs,.DS_Store,annonymize

black:
	@if command -v black > /dev/null; then \
		echo "Running black"; \
		black --check seek; \
	else \
		echo "black not found, please install it!"; \
		exit 1; \
	fi;
	@echo "black passed"

check:
	@$(MAKE) -k black pydocstyle codespell-error

