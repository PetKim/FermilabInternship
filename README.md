# Fermilab Internship for SBND-PRISM

A combination of Python and ROOT scripts I've written.

Frameworks for analayzing different neutrino interaction models utilizing SBND-PRISM.

IMPORTANT: The Python script for the analysis and plots needs heavy optimization as much of the work is straight brute force. There is work that needs to be done to improve the structure and will allow for cleaner analysis in the future. I unfortunately did not have the time nor sufficient knowledge to improve upon it during this time.

For the Python script, it's probably better to stay away from using Pandas for data separation and just stick with NumPy. Other than that, the only thing that needs to be changed are the files being read and what the folder paths are.

The C++ script is similar in the fashion that all you need to do it swap out the ROOT file it reads
