# Documentation

This repository hosts all files for the 6GTandem repository.
All files in `sub-THz-fiber-model` are written in MATLAB and called from the main `simulator` class.


## Installation

<details>
 

 
  <summary>Click me for more details</summary>

 Given that the project consistent of MATLAB and Python, the following should be installed:
- Python
- Matlab
- MATLAB Engine API for Python

### Python installation


### MATLAB installation

Required toolboxes:
- MATLAB
- Signal Processing Toolbox
- Communications Toolbox   
- Phased Array System Toolbox


How to check which toolboxes are required (e.g., to run `do.m`):
```matlab
[fList,pList] = matlab.codetools.requiredFilesAndProducts('do.m')
{pList.Name}'
```


### MATLAB Engine API for Python

[Install MATLAB Engine API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)

Simplified steps:
 - `python -m pip install matlabengine`
 - In Python, import the engine:
```python
import matlab.engine
eng = matlab.engine.start_matlab()
```
</details>

## Models

| **Model category**                       | **Model name**                       | **Included**               | **location** |
|-------------------------------------------------------|-------------------------------------------------------|-------------------------------------------------------|-------------------------------------------------------|
| 6GTandem system simulator                 |  
||    Simulator     | 
||    Modulated signal source and sink     |   
| Stripe model                              |   - | -| `sub-THz-fiber-model`|
|| Polymer microwave/sub-THz fiber           |   
|| Booster units                             |   
|| Radio/Booster unit amplifiers             |   
|| Fiber characteristics                   |   
|| RF switches/splitters/combiners           |   
|| Fiber couplers                            |   
|| Radio units                               |   
|| Phase shifters                            |   
|| Radio unit antennas                       |   
|| Antenna-fiber isolation (cross-talk)      |  
|| Power consumption (PAs and LNAs)          |   
|| CUs                                  		|  
|| Amplifiers                                |   
|| ADC/DAC                                   |   
|| Oscillator                                |   
|| Mixer                                     |   
|| Digital BB processing (incl. source/sink) |   
|| Power consumption                         |   
| Wireless Propagation |
|| Sub-THz wireless channel                  |   
|| Sub-10GHz wireless channel                |   
|| Combined sub-THz/sub-10GHz                |   
| Traffic and mobility    |   
|| Traffic models for different use cases    |   
|| XR Mobility Dataset                       |   
| UEs                       |   
|| Antennas                                  |   
|| Power consumption (UEs)              		|        
