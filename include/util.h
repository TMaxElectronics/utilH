#ifndef UTIL_INC
#define UTIL_INC
    
#include <stdint.h>

#define CONFIG_MAX_LINE_SIZE 128
#define CONFIG_MAX_LINE_COUNT 128

//a handy define to get the size of an array
#define arraySize(X) (sizeof(X) / sizeof(X[0]))

#if !__is_compiling || __has_include("ff.h")
#include "ff.h"

char * CONFIG_getKey(FIL * file, char * keyToFind);
#endif

#define PWL_getRowData(PWL, ROW) (&(PWL->data[ROW * (PWL->preComputedDerivative ? 3 : 2)]))

typedef struct{
    uint32_t listSizeRows; 
    uint32_t preComputedDerivative;
    uint32_t preciceDerivative;
    int32_t * data;
} Pwl_t;

//calculates the datasize for a pwl from values in the descriptor
#define PWL_getDataSize(pwl) (pwl->listSizeRows * sizeof(int32_t) * (pwl->preComputedDerivative ? 3 : 2))

int32_t PWL_getY(int32_t x, Pwl_t * pwl);
void PWL_delete(Pwl_t * pwl, uint32_t freeData);
Pwl_t * PWL_create(int32_t * data, uint32_t rowCount, uint32_t preComputedDerivative, uint32_t preciceDerivative);







typedef struct{
    float R0;
    float T0;
    float Beta;
} NTC_Coefficients_t;

typedef enum{ NTC_MILLI_KELVIN, NTC_MILLI_DEG_CELSIUS, NTC_MILLI_DEG_FAHRENHEIT} NTC_TemperatureUnit_t;

static int32_t NTC_kelvinToUnit(float temperature_K, NTC_TemperatureUnit_t unit);
static float NTC_unitToKelvin(int32_t temperature, NTC_TemperatureUnit_t unit);
int32_t NTC_getTemperatureAtResistance(NTC_Coefficients_t * coefficients, float resistance, NTC_TemperatureUnit_t unit);
float NTC_getResistanceAtTemperature(NTC_Coefficients_t * coefficients, int32_t startTemperature, NTC_TemperatureUnit_t unit);
Pwl_t * NTC_generatePWL(NTC_Coefficients_t * coefficients, int32_t startTemperature, int32_t endTemperature, uint32_t pointCount, NTC_TemperatureUnit_t unit);







int32_t atoiFP(char * a, uint32_t strlen, int32_t baseExponent, uint32_t ignoreUnit);







int32_t qSin(int32_t x);







uint32_t isAsciiNumber(char c);
uint32_t isAsciiSpecialCharacter(char c);








//an array containing bit masks with the correct offset
extern const uint32_t util_bitMasks[];

typedef enum {LITTLE_ENDIAN, BIG_ENDIAN} Endiannes_t;


//gets a 16bit word from an array of bytes
//WARNING: user must make sure that offset+1 bytes are available in the data buffer!
uint16_t get16BitWord(uint8_t * data, uint32_t offset, Endiannes_t endian);
//sets a 16bit word in an array of bytes
//WARNING: user must make sure that the data buffer is large enough!
uint16_t set16BitWord(uint8_t * data, uint32_t offset, Endiannes_t endian, uint32_t value);

//get a bit out of a number of bits in a byte array
//WARNING: user must make sure that enough bits are in the data buffer!
uint32_t getBit(uint8_t * data, uint32_t bitNumber);
//set a bit in a byte array
//WARNING: user must make sure that the data buffer is large enough to even have the bit in question!
void setBit(uint8_t * data, uint32_t bitNumber, uint32_t value);

//divide rounding up the result
//WARNING: x+y must by <= INT32_MAX
int32_t ceil_div(int32_t x, int32_t y);

#endif