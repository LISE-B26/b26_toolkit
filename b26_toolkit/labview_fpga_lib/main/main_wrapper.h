/*
 * Fpga.h
 *
 *  Created on: 08.10.2014
 *      Author: ehebestreit
 *  Edited on 13.01.2016 
 * 		Author JanGieseler
 */
 
#include "NiFpga_main_FPGA.h"
#include <stdio.h>


#ifndef FPGA_H_
#define FPGA_H_

void reset_fpga(NiFpga_Session* session, NiFpga_Status* status);
void start_fpga(NiFpga_Session* session, NiFpga_Status* status);
void stop_fpga(NiFpga_Session* session, NiFpga_Status* status);

// =====================================================================================
// main control parameters
 // =====================================================================================
void set_run_mode(uint16_t value, NiFpga_Session* session, NiFpga_Status* status);
uint16_t read_run_mode(NiFpga_Session* session, NiFpga_Status* status);
void read_stop_all(_Bool state, NiFpga_Session* session, NiFpga_Status* status);
_Bool read_stop_all(NiFpga_Session* session, NiFpga_Status* status);
_Bool read_executing_subvi(NiFpga_Session* session, NiFpga_Status* status);
void set_count_ms(uint32_t value, NiFpga_Session* session, NiFpga_Status* status);
int32_t read_elements_written_to_dma(NiFpga_Session* session, NiFpga_Status* status);

// =====================================================================================
// galvo scan parameters
// =====================================================================================
void set_Nx(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
int16_t read_Nx(NiFpga_Session* session, NiFpga_Status* status);
void set_Vmin_x(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_dVmin_x(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_Ny(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
int16_t read_Ny(NiFpga_Session* session, NiFpga_Status* status);
void set_Vmin_y(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_dVmin_y(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_scanmode_x(uint8_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_scanmode_y(uint8_t state, NiFpga_Session* session, NiFpga_Status* status);
void set_detector_mode(uint8_t state, NiFpga_Session* session, NiFpga_Status* status);
void set_settle_time(uint16_t value, NiFpga_Session* session, NiFpga_Status* status);
uint16_t read_settle_time(NiFpga_Session* session, NiFpga_Status* status);
void set_meas_per_pt(uint16_t value, NiFpga_Session* session, NiFpga_Status* status);
uint16_t read_meas_per_pt(NiFpga_Session* session, NiFpga_Status* status);


// =====================================================================================
// read inputs and outputs
// =====================================================================================
int16_t read_AI0(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI1(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI2(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI3(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI4(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI5(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI6(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AI7(NiFpga_Session* session, NiFpga_Status* status);

int16_t read_AO0(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO1(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO2(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO3(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO4(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO5(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO6(NiFpga_Session* session, NiFpga_Status* status);
int16_t read_AO7(NiFpga_Session* session, NiFpga_Status* status);


void set_AO0(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO1(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO2(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO3(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO4(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO5(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO6(int16_t value, NiFpga_Session* session, NiFpga_Status* status);
void set_AO7(int16_t value, NiFpga_Session* session, NiFpga_Status* status);

_Bool read_DIO0(NiFpga_Session* session, NiFpga_Status* status);
_Bool read_DIO1(NiFpga_Session* session, NiFpga_Status* status);
_Bool read_DIO2(NiFpga_Session* session, NiFpga_Status* status);
_Bool read_DIO3(NiFpga_Session* session, NiFpga_Status* status);

void set_DIO4(_Bool state, NiFpga_Session* session, NiFpga_Status* status);
void set_DIO5(_Bool state, NiFpga_Session* session, NiFpga_Status* status);
void set_DIO6(_Bool state, NiFpga_Session* session, NiFpga_Status* status);
void set_DIO7(_Bool state, NiFpga_Session* session, NiFpga_Status* status);



// =====================================================================================
// ====== FIFO ===
// =====================================================================================
size_t configure_FIFO(size_t requestedDepth, NiFpga_Session* session, NiFpga_Status* status);
void start_FIFO(NiFpga_Session* session, NiFpga_Status* status);
void stop_FIFO(NiFpga_Session* session, NiFpga_Status* status);
void read_FIFO(uint32_t* input, size_t size, NiFpga_Session* session, NiFpga_Status* status,size_t* elementsRemaining);

#endif /* FPGA_H_ */
