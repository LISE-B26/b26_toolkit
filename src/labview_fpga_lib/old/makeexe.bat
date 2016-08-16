@echo off

call gcc -c test_Fpga_c.c NiFpga.c
call gcc -o TEST.exe test_Fpga_c.o NiFpga.o

echo TEST.exe successfully created 

goto:eof


