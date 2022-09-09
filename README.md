# RS-LMTO-ASA

This is a project that has the goal to rewrite the RS-LMTO-ASA code as a OBB modern fortran, using namelists as inputs. This project is inspired by well known electronic structure codes like [Quantum Espresso](https://www.quantum-espresso.org/) and [DyNaMol](https://github.com/araven/DyNaMol). It is a not yet finished project. For the full functioning code, please visit https://gitlab.com/aibman/rslmto.

The RS-LMTO-ASA code is an ab initio code based on Linear Muffin-tin Orbitals within the Atomic Sphere Approximation. 
Instead of solving the eigenvalue problem by diagonalization, RS-LMTO-ASA uses the Haydock recursion method to calculate the local density of states in real-space.


## Installation instructions

1. Download the code 
```bash
   git clone https://github.com/ramoncardias/rs-lmto-asa.git
```

2. Build the code
```bash
   cd rslmto/new_source
   make 
```

## Development instructions

1. Make your changes

2. ***Test*** the code

3. Review the list of changed files 

```
   git status
```
4. Add the changed files to the new revision

```
   git add file1 file2 ...
```
5. Commit the changes

```
   git commit -m "Description of your changes"
```
6. Push the code to `gitlab`

```
   git push
```


## How To Use

The main program is (by default) `rs.x`, it accepts arguments as follows:

First of all, you need to run this executable on a folder that contains input 'namelist-like' file (look at the `example` folder).

1. If you run `./rs.x`, the program will look for `input.nml` to read as input

2. If you want to read a different input file (lets say `my-input.nml`), you can provide the file name as first argument:

```bash
./rs.x my-input.nml
```

3. After run the code for the first time, you will see that the `output` folder was created. The code copy every thing he needs to this folder and run there. If you want it to ran on a different folder, just provide `output=my-output` after the new input name.

```bash
./rs.x input.nml output=my-output
```

4. As you can see in `example` folder, the `input.nml` contain a lot of namelists (json-like structres). If you want to use one of these namelists from another input file, just add in the argument list which one you want and from which file it is.

```bash
./rs.x input.nml atoms=atoms-1.nml
```

Here, the code merge the namelist `atoms` from files `input.nml` and `atoms-1.nml`, last namelists on argument list will have higher priority when common fields emerge. So, if you just wanna add any new field from `atoms-1.nml` but left those from `input.nml` the same, just do:

```bash
./rs.x input.nml atoms=atoms-1.nml atoms=input.nml
```

You can repeat the same namelist as much as you need.

You can do the same with another namelists (like `control`, `lattice`, `self`, `charge`, `calculation`, ...) and mixute them along the argument list, just ensure that:
   - the file contains this namelist
   - there are only one `output=*`
   - there are no spaces in the file names
