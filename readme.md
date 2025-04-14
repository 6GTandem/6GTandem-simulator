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
- Antenna Toolbox


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

| **Model category**   | **Model name**                                    | **Included**      | **location**         | **Python ported** | **Unittested** |
| -------------------- | ------------------------------------------------- | ----------------- | -------------------- | ----------------- | -------------- |
| System-Level         |                                                   |                   |                      |                   |                |
|                      | 6GTandem system simulator (1.5)                   |                   |                      |                   |                |
| Stripe               |                                                   |                   | `sub-THz-link-model` |                   |                |
|                      | Stripe model (2.1.2)                              |                   |                      |                   |
|                      | Polymer microwave/sub-THz fiber (2.2)             | yes               |                      | yes               | no             |
|                      | Booster units (2.3)                               | yes (link)        |                      | yes               | no             |
|                      | Radio/Booster unit amplifiers (2.3.1)             | yes               |                      | yes               | no             |
|                      | RF switches/splitters/combiners (2.3.2)           |                   |                      |                   |                |
|                      | Fiber couplers (2.3.3)                            | yes               |                      | yes               | yes            |
|                      | Radio units (2.4)                                 |                   |                      |                   |                |
|                      | Phase shifters (2.4.1)                            |                   |                      |                   |                |
|                      | Radio unit antennas (2.4.2)                       |                   |                      |                   |                |
|                      | Antenna-fiber isolation (cross-talk) (2.4.3)      |                   |                      |                   |                |
|                      | Power consumption (PAs and LNAs) (2.4.4)          |                   |                      |                   |                |
|                      | Central Unit (CU) (2.5)                           | yes (transmitter) |                      | yes               | no             |
|                      | Amplifiers (2.5.1)                                | yes               |                      | yes               | no             |
|                      | ADC/DAC (2.5.2)                                   | DAC only          |                      | yes (DAC)         | no             |
|                      | Oscillator (2.5.3)                                | yes               |                      | yes               | no             |
|                      | Mixer (2.5.4)                                     | yes               |                      | yes               | yes            |
|                      | Digital BB processing (incl. source/sink) (2.5.5) |                   |                      |                   |                |
|                      | Power consumption (2.5.6)                         |                   |                      |                   |                |
| Wireless Propagation |                                                   |                   |                      |                   |                |
|                      | Sub-THz wireless channel (3.1)                    |                   |                      |                   |                |
|                      | Sub-10GHz wireless channel (3.2)                  |                   |                      |                   |                |
|                      | Combined sub-THz/sub-10GHz (3.3)                  |                   |                      |                   |                |
| Traffic and mobility |                                                   |                   |                      |                   |                |
|                      | Traffic models for different use cases (4.2)      |                   |                      |                   |                |
|                      | XR Mobility Dataset (4.3)                         |                   |                      |                   |                |
| User Equipment (UE)  |                                                   |                   |                      |                   |                |
|                      | Antennas (5.3)                                    |                   |                      | no                | no             |
|                      | Power consumption (UE) (5.4)                      |                   |                      |                   |                |
