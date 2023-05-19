# Raman_Manuscript_2023

## Introduction 
This repository contains the code used to quantify the mechanical behavior of optically stimulated muscle sheets as described in the manuscript "Magnetic Hydrogel Actuators as Mechanically Active Cell Culture Substrates" by Brandon Rios, Angel Bu, Tara Sheehan, Hiba Kobeissi, 
Emma Lejeune, and Ritu Raman (link forthcoming). 

## Code description
The ``code`` folder in this repository contains $4$ python files: ``prepare_files.py``, ``run_code.py``, ``image_analysis.py``, and ``strain_analysis.py``. 

* ``prepare_files.py`` file contains the script required to perform the necessary video pre-processing steps and skip the frames where light stimulation takes place.  
* ``run_code.py`` file runs the tracking pipeline to output and visualize the resutls that mainly include full-field displacements and sub-domain averaged strains. 
* ``image_analysis.py`` script contains the functions to identify fiducial markers, track them across consecutive frames, and compute full-field displacement results. 
* ``strain_analysis.py`` script contains the functions to divide the muscle sheet into a grid box and compute the average Green-Lagrange strain in each of these subdomains.  

For detailed information on the main functionalities contained in ``image_analysis.py`` and ``strain_analysis.py``, we refer the interested user to the main GitHub repository from which this current one is adapted: [MicroBundleCompute](https://github.com/HibaKob/MicroBundleCompute) and the accompanying manuscripty (link forthcoming) 

## Running the code
In the ``tutorial`` folder, we include a single experimental example, "control1.tif", to demonstrate the required initial folder structure as well as the steps to successfully run the provided scripts in the ``code`` folder.

The ``tutorial`` folder should, before running the code, contain the following:

```bash
|___ tutorial
|        |___ original_files
|                |___ "control1.tif"

```


## References to related work 
Related work can be found here:
* Das, S. L., Sutherland, B. P., Lejeune, E., Eyckmans, J., & Chen, C. S. (2022). Mechanical response of cardiac microtissues to acute localized injury. American Journal of Physiology-Heart and Circulatory Physiology, 323(4), H738-H748.

Related repositories include:
* https://github.com/elejeune11/Das-manuscript-2022
* https://github.com/HibaKob/MicroBundleCompute
* https://github.com/HibaKob/SyntheticMicroBundle (synthetic dataset)

## Contact information
For additional information, please contact Emma Lejeune `elejeune@bu.edu` or Hiba Kobeissi `hibakob@bu.edu`.
