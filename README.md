# DLS_GUI

Citation
--------

To cite DLS GUI in publications use
> Salazar, M., Srivastav, H., Srivastava, A., & Srivastava, S. (2023). A User-friendly Graphical User Interface for Dynamic Light Scattering Data Analysis. https://doi.org/10.26434/chemrxiv-2023-v9kt4

Installation
------------

To compile the DLS GUI you will need [Python 3](https://www.python.org/downloads/). 

We also suggest using a virtual environment to keep the compilation clean and isolated. The `venv` package provides an easy way to setup virtual environments and you can find the instructions [here](https://docs.python.org/3/library/venv.html). To install venv package and setup a virtual enviroment named _myenv_:

>`pip install venv`

>`python -m venv myenv`

To activate the virtual environment:
>`myenv\Scripts\activate.bat` (Windows)

>`source myenv/bin/activate` (Unix/Mac)

This package requires `numpy`, `scipy` `matplotlib`, `cz-freeze`, `statsmodels`, and `tkinter` packages. Install these by invoking the `requirements.txt` file provided. Navigate to home directory and execute the following command:

>`pip install -r requirements.txt`

Compile Executable GUI
------------

We are now ready to compile the GUI. In the home directory, execute the following command:
>`python setup.py build_exe --excludes=matplotlib.tests,numpy.random._examples --no-compress` (Windows)

>`python setup.py bdist_mac` (Mac)

There should be new subdirectory _build_ in the home directory containing the executable GUI.  

License
-------

This project is licensed under the [MIT
license](http://en.wikipedia.org/wiki/MIT_License). The documents are
licensed under [Creative Commons Attribution
License](http://creativecommons.org/licenses/by/4.0/).
