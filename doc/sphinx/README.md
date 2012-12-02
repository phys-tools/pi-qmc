This directory has documentation that uses sphinx http://sphinx-doc.org

The language is ReStructureText with some sphinx additions.
Where possible, we link to examples in the system
integration tests (`pi-qmc/doc/system`).

When the code is checked in to github, the http://readthedocs.org
site autogenerates docs at http://pi-qmc.readthedocs.org

You can generate your own local docs with:

    make html

or 

    make latexpdf

The documents will be built in `_build/html/` and `_build/latex/pi-qmc.pdf`.
