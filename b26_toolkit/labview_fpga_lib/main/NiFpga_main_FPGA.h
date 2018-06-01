/*
 * Generated with the FPGA Interface C API Generator 15.0.0
 * for NI-RIO 15.0.0 or later.
 */

#ifndef __NiFpga_main_FPGA_h__
#define __NiFpga_main_FPGA_h__

#ifndef NiFpga_Version
   #define NiFpga_Version 1500
#endif

#include "NiFpga.h"

/**
 * The filename of the FPGA bitfile.
 *
 * This is a #define to allow for string literal concatenation. For example:
 *
 *    static const char* const Bitfile = "C:\\" NiFpga_main_FPGA_Bitfile;
 */
#define NiFpga_main_FPGA_Bitfile "C:\\Users\\Experiment\\PycharmProjects\\b26_toolkit\\src\\labview_fpga_lib\\main\\NiFpga_main_FPGA.lvbitx" 

/**
 * The signature of the FPGA bitfile.
 */
static const char* const NiFpga_main_FPGA_Signature = "A4193EC87477B40F7095B37824BD6B53";

typedef enum
{
   NiFpga_main_FPGA_IndicatorBool_Connector1DIO0 = 0x9E,
   NiFpga_main_FPGA_IndicatorBool_Connector1DIO1 = 0x9A,
   NiFpga_main_FPGA_IndicatorBool_Connector1DIO2 = 0x96,
   NiFpga_main_FPGA_IndicatorBool_Connector1DIO3 = 0x92,
   NiFpga_main_FPGA_IndicatorBool_executingsubscript = 0x12,
   NiFpga_main_FPGA_IndicatorBool_galvo_scan_finished = 0x1A,
} NiFpga_main_FPGA_IndicatorBool;

typedef enum
{
   NiFpga_main_FPGA_IndicatorI16_Connector1AI0 = 0x6E,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI1 = 0x6A,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI2 = 0x66,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI3 = 0x62,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI4 = 0x5E,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI5 = 0x5A,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI6 = 0x56,
   NiFpga_main_FPGA_IndicatorI16_Connector1AI7 = 0x52,
} NiFpga_main_FPGA_IndicatorI16;

typedef enum
{
   NiFpga_main_FPGA_IndicatorI32_datasenttoDMA = 0x1C,
} NiFpga_main_FPGA_IndicatorI32;

typedef enum
{
   NiFpga_main_FPGA_ControlBool_Connector1DIO4 = 0xAE,
   NiFpga_main_FPGA_ControlBool_Connector1DIO5 = 0xAA,
   NiFpga_main_FPGA_ControlBool_Connector1DIO6 = 0xA6,
   NiFpga_main_FPGA_ControlBool_Connector1DIO7 = 0xA2,
   NiFpga_main_FPGA_ControlBool_StopAll = 0xB6,
   NiFpga_main_FPGA_ControlBool_stopmain = 0x16,
} NiFpga_main_FPGA_ControlBool;

typedef enum
{
   NiFpga_main_FPGA_ControlU8_detector_mode = 0x22,
   NiFpga_main_FPGA_ControlU8_scanmodex = 0x26,
   NiFpga_main_FPGA_ControlU8_scanmodey = 0x2A,
} NiFpga_main_FPGA_ControlU8;

typedef enum
{
   NiFpga_main_FPGA_ControlI16_Connector1AO0 = 0x8E,
   NiFpga_main_FPGA_ControlI16_Connector1AO1 = 0x8A,
   NiFpga_main_FPGA_ControlI16_Connector1AO2 = 0x86,
   NiFpga_main_FPGA_ControlI16_Connector1AO3 = 0x82,
   NiFpga_main_FPGA_ControlI16_Connector1AO4 = 0x7E,
   NiFpga_main_FPGA_ControlI16_Connector1AO5 = 0x7A,
   NiFpga_main_FPGA_ControlI16_Connector1AO6 = 0x76,
   NiFpga_main_FPGA_ControlI16_Connector1AO7 = 0x72,
   NiFpga_main_FPGA_ControlI16_N_x = 0x36,
   NiFpga_main_FPGA_ControlI16_N_y = 0x42,
   NiFpga_main_FPGA_ControlI16_Vmin_x = 0x3E,
   NiFpga_main_FPGA_ControlI16_Vmin_y = 0x4A,
   NiFpga_main_FPGA_ControlI16_dVmin_x = 0x3A,
   NiFpga_main_FPGA_ControlI16_dVmin_y = 0x46,
} NiFpga_main_FPGA_ControlI16;

typedef enum
{
   NiFpga_main_FPGA_ControlU16_measurements_per_pt = 0x2E,
   NiFpga_main_FPGA_ControlU16_run_mode = 0xB2,
   NiFpga_main_FPGA_ControlU16_settle_time_us = 0x32,
} NiFpga_main_FPGA_ControlU16;

typedef enum
{
   NiFpga_main_FPGA_ControlU32_CountmSec = 0x4C,
} NiFpga_main_FPGA_ControlU32;

typedef enum
{
   NiFpga_main_FPGA_TargetToHostFifoI32_DMA = 0,
} NiFpga_main_FPGA_TargetToHostFifoI32;

#endif
