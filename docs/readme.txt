This folder contains the source documents for the github pages.

The github pages lives in the orphaned branch gh-pages.

To build, commands like these are needed:

``` shell
deactivate # in a virtual environment is in use
python -m venv .venv
source .venv/bin/activate
python -m pip install sphinx sphinx_rtd_theme
make html
```

The pages were intialized with

``` shell
sphinx-quickstart
```
