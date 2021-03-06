# pcNN_solid
Add-on to install pcNN XC functional into VASP source code. 

This add-on is available only to owners of a valid VASP license.


## requirements
VASP 5.4.4 with source code (It does not work with VASP 6.)  
Python3  
PyArmor  


## Installation
1. Backup "vasp/src/metagga.F" to another directory.
2. Run "patch_implement_pcNN.py". You will be asked the absolute path to metagga.F.
3. "metagga_patched.F" will be generated. Rename it to "metagga.F" and replace "vasp/src/metagga.F" with it.
4. Compile VASP.

## Usage 
1. Place "nnparams" folder into the same directory as INCAR.

(example)
```
(target folder) -- nnparams
                -- INCAR
                -- POSCAR
                -- KPOINTS
                -- POTCAR
                ...
```

2. Add (or modify) keyword "METAGGA = NNMGGA" to INCAR.
3. Run your calculation with the compiled VASP.


## Citing this XC functinoal
@misc{nagai2021machinelearningbased,
      title={Machine-Learning-Based Exchange-Correlation Functional with Physical Asymptotic Constraints}, 
      author={Ryo Nagai and Ryosuke Akashi and Osamu Sugino},
      year={2021},
      eprint={2111.15593},
      archivePrefix={arXiv},
      primaryClass={cond-mat.mtrl-sci}
}


## Notice
This add-on is distributed under an agreement from VASP Software GmbH.  

- This add-on comes without any waranty.  
- If the add-on is modified by inclusion of additional source code lines from the main VASP source code, 
the developer must send the revised version of the add-on to VASP Software GmbH for further review.


## License
Copyright 2021 Ryo Nagai

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

