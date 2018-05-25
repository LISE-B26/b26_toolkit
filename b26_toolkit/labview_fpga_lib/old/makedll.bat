cls
@echo off
REM this batch file compiles the .c file (passed as parameter) and creates a .dll

REM if [%1]==[] goto usage
REM if NOT %~x1==.c goto wrongfiletype

REM call gcc -c %1 NiFpga.c
REM call gcc -shared -o %~n1.dll %~n1.o

REM ==== complile for simple inputs and outputs =====
REM call gcc -c Fpga.c NiFpga.c
REM call gcc -shared -o FPGAlib.dll Fpga.o NiFpga.o
REM echo FPGAlib.dll successfully created 

REM ==== complile for PID =====
call gcc -c FPGA_PID_Wrapper.c NiFpga.c
call gcc -shared -o FPGA_PID_lib.dll FPGA_PID_Wrapper.o NiFpga.o
echo FPGAlib.dll successfully created 



goto:eof

:wrongfiletype
@echo this is not a .c file
:usage
@echo please provide a source .c file
@echo usage: %0 source.c

