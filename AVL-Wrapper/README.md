# AVL-Wrapper
### Aero team AVL wrapper
ADaX stands for AVL Designer and Executor. It is a Python library and mini-interface for creating ".avl" and ".run" files and automatically running them through AVL.

### Structure
"ADaX.py" is the main script to run. It has three options: Reload_Geometry, Run_AVL, and Autodebug. 
- Reload_Geometry: The geometry is handled in the file "Code/\_\_Aircraft_Geometry__.py", where surfaces, sections, and controls may be fully defined using Python variables. If Reload_Geomery is True, then the aircraft geometry script will be run before running AVL. The script allows a unique file path to be created when writing the geometry, but reloading the geometry is optional to avoid overwriting the current geometry.
- Run_AVL: Chooses whether or not the custom run case is executed. The custom run case is set up in "Code/\_\_AVL_Commands__.py". This may be turned off to purely create the aircraft geometry without running AVL.
- Autodebug: Set to True for an experimental autodebugger to notify the user of common errors. If output results do not reflect the run case you are trying to execute, there is likely an error preventing flow execution or failure to trim.

"CaseMaker.py" is the file for creating run cases. It reflects the layout of AVL's ".run" files and is capable of creating new ".run" files, append to existing ones, or overwriting existing cases inside a ".run" file

### Documentation
Most documentation is included in the file "Code/Codebase.py", explaining commands used in "Code/\_\_AVL_Commands__.py". However, documentation is a WIP. I've gotten scarce feedback from sharing the code with others in senior design, and I've made some changes based on how I saw them use it. Therefore, I'd like some fresh outside perspectives on the code's layout and what is still required for documentation (@Esther @Nick!!!)

### Getting Started
Open AVL-Wrapper-main in Visual Studio, as all the code's directories are based on this folder being the working directory. Maybe someone in the future could do some directory management. 
