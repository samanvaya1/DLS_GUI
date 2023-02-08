# DLS-GUI

Background
--------

[Dynamic light scattering](https://en.wikipedia.org/wiki/Dynamic_light_scattering) (DLS) is a commonly used analytical tool for characterizing the size distribution of colloids in a dispersion or a solution. Typically, the intensity of a scattering produced from the sample at a fixed angle from an incident laser beam is recorded as a function of time and converted into time autocorrelation data, which can be inverted to estimate the distribution of colloid diffusivity to estimate the colloid size distribution.

Introduction
--------

DLS-GUI is a user-friendly graphical user interface (GUI) for analyzing the measured scattering intensity time autocorrelation data using both the cumulant expansion method and regularization methods, with the latter implemented using various commonly employed algorithms including NNLS, CONTIN, REPES, and DYNALS. Additionally, the GUI also enables a comparison of the size distributions generated by various algorithms and an evaluation of their performance. To learn more about the software and/or download the compiled version please visit the [Lab Website](https://sites.google.com/view/srivastava-lab/software).

Authors: Matthew Salazar, Harsh Srivastav, Abhishek Srivastava, Samanvaya Srivastava  
[Self-Assembly of Soft Materials Laboratory](https://sites.google.com/view/srivastava-lab), University of California, Los Angeles

To cite DLS GUI in publications please use:
> Salazar M, Srivastav H, Srivastava A, Srivastava S. A User-friendly Graphical User Interface for Dynamic Light Scattering Data Analysis. ChemRxiv. Cambridge: Cambridge Open Engage; 2023 https://doi.org/10.26434/chemrxiv-2023-v9kt4  

DLS GUI Preview:  
<img src="https://user-images.githubusercontent.com/24664516/217431803-7c9c7a38-d4ad-4ffe-b928-2503e9d097bb.png" alt="GUI_Mac" width=50% height=50%/>

Installation
------------

To compile the DLS GUI you will need to install [Python 3](https://www.python.org/downloads/). 

We also suggest using a virtual environment to keep the compilation clean and isolated. The `venv` package provides an easy way to setup virtual environments and you can find the detailed description and instructions [here](https://docs.python.org/3/library/venv.html).

To install venv package and setup a virtual enviroment named _myenv_:

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

There should be a new subdirectory _build_ in the home directory containing the executable GUI.  

Using the GUI
-------
Detailed descriptions of the features along with instructions for data loading, parameter adjustment, analysis, algorithms comparison, and exporting are provided in the attached manual _DLS GUI Manual.pdf_.  
<img src="https://user-images.githubusercontent.com/24664516/217431770-9b091b5b-6803-4afb-90e6-5af40e07dbb4.JPG" alt="Manual_Contents" width=40% height=40%/>

License
-------

This project is licensed under the [MIT
license](http://en.wikipedia.org/wiki/MIT_License). The documents are
licensed under [Creative Commons Attribution
License](http://creativecommons.org/licenses/by/4.0/).

References
-------
1. Salazar M, Srivastav H, Srivastava A, Srivastava S. A User-friendly Graphical User Interface for Dynamic Light Scattering Data Analysis. ChemRxiv. Cambridge: Cambridge Open Engage; 2023 [[Link](https://doi.org/10.26434/chemrxiv-2023-v9kt4)]
2. Scotti AE, Liu W, Hyatt JS, Herman ES, Choi HS, Kim JW, Lyon LA, Gasser U, Fernandez-Nieves A. The CONTIN algorithm and its application to determine the size distribution of microgel suspensions. The Journal of chemical physics. 2015 Jun 21;142(23):234905 [[Link](https://aip.scitation.org/doi/abs/10.1063/1.4921686)]
3. Berne BJ, Pecora R. Dynamic light scattering: with applications to chemistry, biology, and physics. Courier Corporation; 2000. [[Link](https://books.google.com/books?hl=en&lr=&id=vBB54ABhmuEC&oi=fnd&pg=PA1&dq=dynamic+light+scattering&ots=L7pGK2qksd&sig=Sp3ud4IpPO4jt1UrWoY8wqNpy_w)]
