# Fermilab Internship: Researching SBND-PRISM

A combination of Python and ROOT scripts I've written (with some help of course).

Frameworks for analayzing different neutrino interaction models utilizing SBND-PRISM.

IMPORTANT: The Python code for the analysis and plots needs heavy optimization as much of the work is straight up brute force. There is work that needs to be done to improve the structure and will allow for cleaner analysis in the future. I unfortunately did not have the time nor sufficient knowledge to improve upon it during this time.

For the Python code, it's probably better to stay away from using Pandas for data separation and just stick with NumPy. Other than that, the only thing that needs to be changed are the files being read and what the folder paths are. Please disregard the error pages you see on the Jupyterlab script. I had to reset my laptop at one point near the end of the internship and it messed with some of the file paths it reads and can be easily fixed.

The C++ script is similar in the fashion that all you need to do it swap out the ROOT file it reads.

It should be noted as well the only difference between the 2 C++ scripts is that a boolean variable has been added in the CCQE file which only extracts events that have been tagged with that specific requirement. It can be swapped out for different interaction types if needed.

I have left all comments in as they serve as a reminder and were very helpful for me during this period.

I will be working on improving it when I have the time.
