<h1 align='center'>sharpInterFoam</h1>


Customized OpenFOAM solver for two-phase flow based on a coupling of the THINC/SW and the S-CLSVOF methods.
The code for the paper: [*Int. J. Microgravity Sci. Appl.*, 38(3) 380301 (2021)](https://doi.org/10.15011/jasma.38.380301)


#### Efficient Implementation of Two-Phase Flow Solver Based on THINC/SW and S-CLSVOF on Unstructured Meshes<br>
Suguru Shiratori<sup>1,*</sup>, 
Takuro Usui<sup>2</sup>,
Shiho Koyama<sup>2</sup>,
Shumpei Ozawa<sup>3</sup>,
Hideaki Nagano<sup>1</sup>,
Kenjiro Shimano<sup>1</sup>
<br>
<sup>1</sup> Dept. Mechanical Systems Engineering, Tokyo City University, Japan <br>
<sup>2</sup> Graduate School of Integrative Science and Engineering, Tokyo City University, Japan <br>
<sup>3</sup> Department of Advanced Materials Science and Engineering, Chiba Institute of Technology, Chiba, Japan <

## Tested OpenFOAM versions
* OpenFOAM-v2412

## Build guide
* Install OpenFOAM-v2412
* Build sharpInterFoam
```bash
$ source $WM_PROJECT_DIR/etc/bashrc
$ cd lib/THINC
$ wmake
$ cd ../SCLSVOF
$ wmake
$ cd ../../solver
$ wmake
```

## Usage
Example cases are in `test_case` directory.

Parameters particular for <b>sharpInterFoam</b> are in the file: `system/fvSolutions`.
```
"alpha.water.*"
{
  SCLSVOF  // parameters for SCLSVOF
  {
    correctPsi      yes;            // Yes: SCLSVOF is enabled, No: curvature is calculated same as interFoam
    densityScaled   yes;            // Switch for density-scaled CSF model
    initializeAtanh yes;            // Method for re-initialization of of level-set function. Yes: atan, No: linear
    densityFunctionHeaviside  yes;  // How to evaluate mixture density. Yes: Heaviside function, No: alpha
    widthFactor     1.5;            // 
    denomDeltaTau   20.0;           //
    factorNumLoop   3.0;            //
    tolMagGradPsi   1e-8;           //
  }
 
  THINC    // parameters for THINC
  {
    enable          yes;            // Yes: THINC is enabled, No: calculated by MULES (interFoam)
    epsAlpha        1e-13;          // 
    epsBeta         1e-3;           //
  }
}
```


## Citation
```bibtex
@article{Shiratori_IJMSA_2021,
  title={Efficient Implementation of Two-Phase Flow Solver Based on THINC/SW and S-CLSVOF on Unstructured Meshes},
  author={Suguru Shiratori and Takuro Usui and Shiho Koyama and Shumpei Ozawa and Hideaki Nagano ad Kenjiro Shimano},
  journal={Int. J. Microgravity Sci. Appl.},
  volume={38},
  number={3},
  pages={380301},
  year={2021},
  doi={10.15011/jasma.38.380301},
  }
```
