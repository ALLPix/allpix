/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *
 */

#define MAX_FRAME_ROW 2048
#define MAX_FRAME_COL 2048
#define MAX_FRAME_OCC 65536

#define META_DATA_LINE_SIZE 512
#define ACTUAL_VALUE_OFFSET 2

#define __4thbyte 3
#define __3thbyte 2
#define __2thbyte 1
#define __1stbyte 0

#define ACQ_MODE_STRING                   "Acq mode"
#define ACQ_MODE_EXTRA                    "(\"Acquisition mode\")"
#define ACQ_MODE_CODE                     10

#define ACQ_TIME_STRING                   "Acq time"
#define ACQ_TIME_EXTRA                    "(\"Acquisition time [s]\")"
#define ACQ_TIME_CODE                     20

#define APPLIED_FILTERS_STRING            "Applied filters"
#define APPLIED_FILTERS_EXTRA             "(\"Names of filters which were applied to frame\")"
#define APPLIED_FILTERS_CODE              30

#define AUTO_ERASE_INTERVAL_STRING        "Auto Erase interval"
#define AUTO_ERASE_INTERVAL_EXTRA         "(\"In triggered mode autoerase matrix after [s]\")"
#define AUTO_ERASE_INTERVAL_CODE          40

#define AUTOERASE_INTERVAL_COUNTER_STRING "Autoerase interval counter"
#define AUTOERASE_INTERVAL_COUNTER_EXTRA  "(\"Autoerase interval counter\")"
#define AUTOERASE_INTERVAL_COUNTER_CODE   50

#define BS_ACTIVE_STRING                  "BS active"
#define BS_ACTIVE_EXTRA                   "(\"Back side preamp. enabled\")"
#define BS_ACTIVE_CODE                    60

#define CHIPBOARDID_STRING                "ChipboardID"
#define CHIPBOARDID_EXTRA                 "(\"Medipix or chipboard ID\")"
#define CHIPBOARDID_CODE                  70

#define COINC_LIVE_TIME_STRING            "Coinc live time"
#define COINC_LIVE_TIME_EXTRA             "(\"Last coincidence live time [s]\")"
#define COINC_LIVE_TIME_CODE              80

#define COINCIDENCE_DELAY_STRING          "Coincidence delay"
#define COINCIDENCE_DELAY_EXTRA           "(\"Maximum coincidence delay (1=min, 256=max [=64us])\")"
#define COINCIDENCE_DELAY_CODE            90

#define COINCIDENCE_MODE_STRING           "Coincidence mode"
#define COINCIDENCE_MODE_EXTRA            "(\"Coinc Mode: On/Off=0x01, TrgInLevel=0x02, Notify=0x04, NotifInitLevel=0x08, NotifyValid=0x10, NotifyByPulse=0x20\")"
#define COINCIDENCE_MODE_CODE             100

#define COUNTERS_STRING                   "Counters"
#define COUNTERS_EXTRA                    "(\"Counts: Triggers, Valid coinc, Invalid coinc\")"
#define COUNTERS_CODE                     110

#define DACS_STRING                       "DACs"
#define DACS_EXTRA                        "(\"DACs values of all chips\")"
#define DACS_CODE                         120

#define HV_STRING                         "HV"
#define HV_EXTRA                          "(\"Bias Voltage [V]\")"
#define HV_CODE                           130

#define HW_TIMER_STRING                   "Hw timer"
#define HW_TIMER_EXTRA                    "(\"Hw timer mode\")"
#define HW_TIMER_CODE                     140

#define INTERFACE_STRING                  "Interface"
#define INTERFACE_EXTRA                   "(\"Medipix interface\")"
#define INTERFACE_CODE                    160

#define MPX_CLOCK_STRING                  "Mpx clock"
#define MPX_CLOCK_EXTRA                   "(\"Medipix clock [MHz]\")"
#define MPX_CLOCK_CODE                    180

#define MPX_TYPE_STRING                   "Mpx type"
#define MPX_TYPE_EXTRA                    "(\"Medipix type (1-2.1, 2-MXR, 3-TPX)\")"
#define MPX_TYPE_CODE                     200

#define POLARITY_STRING                   "Polarity"
#define POLARITY_EXTRA                    "(\"Detector polarity (0 negative, 1 positive)\")"
#define POLARITY_CODE                     220

#define START_TIME_STRING                 "Start time"
#define START_TIME_EXTRA                  "(\"Acquisition start time\")"
#define START_TIME_CODE                   240

#define START_TIME_S_STRING               "Start time (string)"
#define START_TIME_S_EXTRA                "(\"Acquisition start time (string)\")"
#define START_TIME_S_CODE                 260

#define TIMEPIX_CLOCK_STRING              "Timepix clock"
#define TIMEPIX_CLOCK_EXTRA               "(\"Timepix clock (0-3: 10MHz, 20MHz, 40MHz, 80MHz)\")"
#define TIMEPIX_CLOCK_CODE                280

#define TRIGGER_TIME_STRING               "Triger time"
#define TRIGGER_TIME_EXTRA                "(\"Last trigger time [s]\")"
#define TRIGGER_TIME_CODE                 300


#define TYPE_256x256_STRING  "256x256_FORMAT"
#define TYPE_XYC_STRING      "XYC_FORMAT"
#define TYPE_XC_STRING       "XC_FORMAT"
#define TYPE_XYC_BIN_STRING  "XYC_BINARY_FORMAT"

#define CHARS_INLINE_256x256_FORMAT  516
#define CHARS_INLINE_XYC_FORMAT      18+5 /* trying to meet the
					     maximun lenght in the
					     case of XYC format. 5 is
					     just an offset */

#define FRAME_TYPE_STRING        "format = "
#define BINARY_FRAME_TYPE_STRING "B000000001"

// The following definitions are used to comply with a particular detector technology
//  and be able to compare with real data.  We need something generic for other cases.
// Format flags for saving to file (MgrSaveFrameType or flags for perform acq. functions)
// flags FSAVE_BINARY/FSAVE_ASCII and flags FSAVE_I16/FSAVE_U32/FSAVE_DOUBLE are mutually exclusive
// if FSAVE_SPARSEX and FSAVE_SPARSEXY is not specified full matrix is saved
#define FSAVE_BINARY        0x0001      // save in binary format
#define FSAVE_ASCII         0x0002      // save in ASCII format
#define FSAVE_APPEND        0x0004      // append data file to existing file if exists
#define FSAVE_I16           0x0010      // save as 16bit integer
#define FSAVE_U32           0x0020      // save as unsigned 32bit integer
#define FSAVE_DOUBLE        0x0040      // save as double
#define FSAVE_NODESCFILE    0x0100      // do not save description file
#define FSAVE_SPARSEXY      0x1000      // save only nonzero position in [x y count] format
#define FSAVE_SPARSEX       0x2000      // save only nonzero position in [x count] format
#define FSAVE_NOFILE        0x8000      // frame will not be saved :)
