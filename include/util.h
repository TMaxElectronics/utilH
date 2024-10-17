#ifndef UTIL_INC
#define UTIL_INC
    
#include <stdint.h>

#define CONFIG_MAX_LINE_SIZE 128
#define CONFIG_MAX_LINE_COUNT 128

#if __has_include("ff.h")
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

#endif