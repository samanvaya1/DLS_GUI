import os
import cx_Freeze
from cx_Freeze import setup, Executable
import sys

packages = ['os', 'matplotlib', 'scipy', 'numpy']
executables = [cx_Freeze.Executable('DLS_GUI.py')]
options = {'build_exe': {"includes": ['numpy', 'tkinter'],
                         'packages': packages,

                         }, }

cx_Freeze.setup(
    name='Test1',
    options=options,
    description='Extraction of data',
    executables=executables
)

# To run must use this line in command window or terminal depending on type of operating system
# python setup.py build_exe --excludes=matplotlib.tests,numpy.random._examples --no-compress
#python setup.py bdist_mac
#Auto Pyinstaller
#CMD: auto-py-to-exe
