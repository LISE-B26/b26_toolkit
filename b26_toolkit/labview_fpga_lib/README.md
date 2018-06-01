Author: Jan Gieseler
Date: Feb 29th 2016


# file structure
- for each labview programm create a folder, e.g. read_ai_ao
- this folder contains
    .the files created by the labview C API creator (FPGA bitfile (.lvbitx) and .c, .h files))
    .the python wrapper file named as the folder, e.g. read_ai_ao_wrapper.py
    .the .c and .h wrapper file names as the folder with appendix _wrapper.c and _wrapper.h, e.g. read_ai_ao_wrapper.c and read_ai_ao_wrapper.h
(June 2nd: following is probably not true)
(- in the main folder (labview_fpga_lib) there is a .py file with the same name as the subfolder, e.g. read_ai_ao.py) 


# complile dll (see also how to compile under windows below!)
1.) compile the labview FGPGA (creates .lvbitx file) 
2.) convert into .h with C API generator (right click in FPGA .vi in the project explorer)
    .creates NiFpga.c, NiFpga.h, NiFpga_PID.c and copies NiFpga_PID.lvbitx to target folder
    .choose subfolder as target (e.g. src/labview_fpga_lib/read_ai_ao)
    .in NiFpga_XXX.h change "#define NiFpga_XXX.lvbitx" to "#define PATH_TO_SRC/labview_fpga_lib/read_ai_ao/NiFpga_XXX.lvbitx" where XXX is the name of the labview file
3.) create/ modify wrapper .c file (e.g. read_ai_ao_wrapper.c and read_ai_ao_wrapper.h)
    .watch out for varialbe names in .c file which contain the folder name such as NiFpga_FPGA_read_ai_ao_IndicatorI16_Connector1AI0
4.) run batch file makedll.bat *folder_name*, e.g. makedll.bat read_ai_ao


# write python wrapper file
1.) import .dll, e.g. _libfpga = WinDLL('C:/Users/Experiment/PycharmProjects/PythonLab/src/labview_fpga_lib/read_ai_ao/read_ai_ao.dll')
2.) access c-functions as following example _libfpga.set_DIO7(value, byref(session), byref(status))


### how to compile under windows ###
install http://mingw-w64.org/ (if for 64 bit windows otherwise go for the 32bit version)
note that python has to be 64bit, too otherwise it will give a windows error that dll is not a valid 32bit file
browse to the folder where you installed minwg-w64 (e.g. C:\mingw-w64\x86_64-6.1.0-posix-seh-rt_v5-rev0)and execute mingw-w64.bat to open command window with correctly set path variables
browse to folder where the makedll.bat file is located (e.g. C:\Users\Experiment\PycharmProjects\b26_toolkit\src\instruments\labview_fpga_lib) and run .bat file








