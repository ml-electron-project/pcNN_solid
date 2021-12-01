# pcNN_solid
VASP version of pcNN XC functional. This add-on is only available to owners of a valid VASP license.


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

2. Add (or modify) the keyword "METAGGA = NNMGGA" to the INCAR.
3. Run the calculation with the VASP program compiled in the Installation section.

## License
This add-on is distributed under an agreement from VASP Software GmbH.  

- This add-on comes without any waranty.  
- Commercial use and patent use are prohibited.
- If the add-on is modified by inclusion of additional source code lines of the Add-on to the VASP software code, 
the developer must send the revised version of the add-on to VASP Software GmbH for further review.


## Citing this XC functinoal
@misc{nagai2021machinelearningbased,
      title={Machine-Learning-Based Exchange-Correlation Functional with Physical Asymptotic Constraints}, 
      author={Ryo Nagai and Ryosuke Akashi and Osamu Sugino},
      year={2021},
      eprint={2111.15593},
      archivePrefix={arXiv},
      primaryClass={cond-mat.mtrl-sci}
}
