#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#if !__is_compiling || __has_include("ff.h")
#include "ff.h"
#endif

#if !__is_compiling || __has_include("TTerm.h")
#include "TTerm.h"
#endif

#include "include/util.h"
#include "System.h"

/*
 * peicewise linear function algorithm, allows for fast lut implementations
 * 
 * usage: call with any x value and information to y LUT you created to get the linear approximation of the graph inbetween the two nearest matching points
 *      
 * preComputedDerivative tells the interpreter wether you have a third dataset in each row containing the derivative at that point or not    
 *      
 * PWL row-format:
 *   preComputedDerivative=0:
 *      {xValue,yValue}
 *      
 *      - xValue: x coordinate of the point
 *      - yValue: y coordinate of the point
 *   preComputedDerivative=1:
 *      {xValue,yValue,dY/dX}
 *      - x/yValue: as in preComputedDerivative=0
 *      - dY/dX: pre-computed rate of change in between the point and the next one in the list
 *  
 *      NOTE: the list must be sorted by x values in ascending order (x[0] < x[1] < x[2]...)
 *      NOTE: if you need this to be fast, make sure to pre-calculate the derivatives as this saves a division for every conversion
 */
int32_t PWL_getY(int32_t x, Pwl_t * pwl){
    //check if we even got a PWL
    if(pwl == NULL) return 0;
    
    if(pwl->listSizeRows < 2){
        //can't approximate any function if all we have is a single point. We need at least two
        //configASSERT(0);
        return 0;
    }
    
    int32_t * currentRow = NULL;
    int32_t * lastRow = NULL;
    
    //find the neighbouring points
    int32_t entry = 0;
    for(entry = 0; entry < pwl->listSizeRows; entry++){
        //calculate data pointer to current row. If preComputedDiff is on then there is a third value in the row (dX/dY between the current and next entry) otherwise only 2
        currentRow = (pwl->data) + entry * (pwl->preComputedDerivative ? 3 : 2);
        
        //check if by any chance x exactly matches the value in the row
        if(currentRow[0] == x){
            //oh yes actually it does => just return the y value
            return currentRow[1];
        }
        
        //check if the entry we are looking at has a x value larger than the one we are looking for
        if(currentRow[0] > x){
            //entry found :)
            
            //check if we are at the first entry in the list
            if(lastRow == NULL){ 
                //yes we are, that means that the x value we are looking for lies outside the boundary of the pwl. Try to interpolate with the first two points
                
                //is there at least one point ahead of us?
                if(entry + 1 >= pwl->listSizeRows){
                    //hmm no, this shouldn't technically be possible as we already check for a list size of at least 2. just return, we can't do any calculations now anyway
                    
                    //this btw means that we are left of the first point in the graph, but no other point exists after the first one (aka we only have one in total)
                    return 0;
                }
                
                //change the pointers to point to the correct points
                lastRow     = currentRow;
                currentRow  = (pwl->data) + (entry + 1) * (pwl->preComputedDerivative ? 3 : 2);
                
                break;
            }
            
            //point found in list, no change to pointers needed
            break;
        }
        
        //entry not yet found, update last point
        lastRow = currentRow;
    }
    
    //check if there was no entry in the list larger than the x we are looking for. If so just use the maximum possible value for calculation.
    if(lastRow == currentRow){
        //this means loop exited due to entry == pwl.listSizeRows. Move lastRow back one entry
        
        //can we actually do that?
        if(entry - 2 < 0){
            //hmm no, this again means there is only one point in the list... can't do anything
            return 0;
        }
        
        lastRow = (pwl->data) + (entry-2) * (pwl->preComputedDerivative ? 3 : 2);
    }
    
    //now interpolate between the two points pointed to by lastRow (start point, left of X) and currentRow (end point, right of x))
            
    //first calculate the x offset from the start point
    int32_t localX = x - lastRow[0];
    int32_t localY = lastRow[1];        

    //calculate derivative of pwl between start and end points
    int32_t dYdX = 0;
    if(pwl->preComputedDerivative){
        //nothing to compute, just read the value from the start point dataset
        dYdX = lastRow[2];
    }else{
        int32_t dy = (currentRow[1] - lastRow[1]) * (pwl->preciceDerivative ? 256 : 1);
        int32_t dx = currentRow[0] - lastRow[0];
        if(dx != 0) dYdX = dy / dx; else dYdX = 0;
    }

    //calculate linear function as 
    //   y = m     *   x       + b         with m=dYdX, x=localX, b=localY
    // option "preciceDerivative" means the derivative value given is multiplied by 256
    return  (pwl->preciceDerivative ? ((dYdX   *   localX) >> 8)  :  (dYdX  *   localX))  + localY;
}

/* 
 * Function to allocate PWL memory
 * 
 * note: 
 *      if data=NULL will cause the function to try to allocate suitably sized memory for it aswell.
 *      You can specify a data pointer if you want something like a dynamically allocated header but a const dataset
 */
Pwl_t * PWL_create(int32_t * data, uint32_t rowCount, uint32_t preComputedDerivative, uint32_t preciceDerivative){
    
//try to allocate pwl header memory
    Pwl_t * pwl = pvPortMalloc(sizeof(Pwl_t));
    if(pwl == NULL) return NULL;    //didn't work => don't even try to continue
    
    pwl->listSizeRows = rowCount;
    pwl->preComputedDerivative = preComputedDerivative;
    pwl->preciceDerivative = preciceDerivative;
    
//do we need to allocate data memory?
    if(data == NULL){
        //yes! => try to do so (*3 is for the row size of three int32's)
        pwl->data = pvPortMalloc(sizeof(int32_t) * (preComputedDerivative ? 3 : 2) * rowCount);
        
        if(pwl->data == NULL){ 
            vPortFree(pwl);
            return NULL;    //didn't work => don't even try to continue
        }
    }else{
        //no. Just assign the data the caller gave us
        pwl->data = data;
    }
    
    return pwl;
}

/* 
 * Function to free a PWLs memory
 */
void PWL_delete(Pwl_t * pwl, uint32_t freeData){
    //is there anything to free?
    if(pwl == NULL) return;
    
    //do we need to free the data? If yes then do so
    if(freeData) vPortFree(pwl->data);
    
    //free the header
    vPortFree(pwl);
}








/*
 * NTC Tool - generates a PWL with a specified number of rows that allows for faster conversion of resistance to temperature
 *      TODO: not sure if we actually need this, but it might be useful: NOTE: for conversion of divider ratio to resistance check NTC_calculateNTCResistance()
 * 
 * usage: 
 *      specify NTC type with the coefficients (R0, T0 and Beta are all needed)
 *      List size is determined by the start & endTemperature (I recommend only specifying the temperature region of interest to you here to increase accuracy) and the point Count. Use at least (TODO: figure out how many points are recommended lol)
 *      select the unit of the PWL with the NTC_TemperatureUnit_t field
 * 
 *      To perform the conversion pass the generated PWL to the PWL_getY() function together with the resistance in Ohms
 */
Pwl_t * NTC_generatePWL(NTC_Coefficients_t * coefficients, int32_t startTemperature, int32_t endTemperature, uint32_t pointCount, NTC_TemperatureUnit_t unit){
//are the parameters valid?
    if(startTemperature >= endTemperature || pointCount < 2 || coefficients == NULL) return NULL;
    
//create a new PWL. Just return if that doesn't work
    Pwl_t * pwl = PWL_create(NULL, pointCount, 1, 1);
    if(pwl == NULL) return NULL;
    
//calculate start and end resistance values. The end value is the one matching the start temperature as we need to sort by ascending resistance
    
    float startResistance = NTC_getResistanceAtTemperature(coefficients, endTemperature, unit);
    float endResistance = NTC_getResistanceAtTemperature(coefficients, startTemperature, unit);
    float resistanceStep = ((float) endResistance - (float) startResistance) / (float) pointCount;
    
    //step through each list entry and calculate the temperature corresponding to the resistance
    //for reference: PWL format with preComputedDerivative=1 is {xValue,yValue,dY/dX}
    
    float currentResistance = startResistance;
    int32_t * lastRow = NULL;
    
    for(int32_t i = 0; i < pointCount; i++){
        //generate a pointer to the current row
        int32_t * currentRow = &(pwl->data[i * 3]);
        
        //first row entry is the x value (in our case the resistance)
        currentRow[0] = (int32_t) currentResistance;
        
        //second row entry is the resulting y value (aka the temperature, converted into the desired unit)
        currentRow[1] = NTC_getTemperatureAtResistance(coefficients, currentResistance, unit);
        
        //third row entry is the derivative between this point and the last
        //this is always calculated for the last entry instead of the current one, so we skip this for the first entry
        //we don't need to add a derivative for the very last row of the PWL as it is never used
        if(i != 0){
            //derivativce is dy/dx
            int32_t dy = (currentRow[1] - lastRow[1]) * 256;
            int32_t dx = currentRow[0] - lastRow[0];
            if(dx != 0) lastRow[2] = dy / dx; else lastRow[2] = 0;
        }
        
        lastRow = currentRow;
        currentResistance += resistanceStep;
    }
    
    return pwl;
}

float NTC_getResistanceAtTemperature(NTC_Coefficients_t * coefficients, int32_t temperature, NTC_TemperatureUnit_t unit){
    //convert temperature to Kelvin
    float t1_K = NTC_unitToKelvin(temperature, unit);
    
    //calculate Resistance according to Beta formula
    // R1 = R0 * exp(Beta * (1/T1 - 1/T0))
    return coefficients->R0 * exp(coefficients->Beta * ((1/t1_K) - (1/coefficients->T0)));
}

int32_t NTC_getTemperatureAtResistance(NTC_Coefficients_t * coefficients, float resistance, NTC_TemperatureUnit_t unit){
    //use the Beta formula to calculate the temperature
    float temperature_K = 1.0 / ( (log(resistance / coefficients->R0) / coefficients->Beta ) + 1.0 / coefficients->T0);
    
    //return the value after converting to the desired unit
    return NTC_kelvinToUnit(temperature_K, unit);
}

static int32_t NTC_kelvinToUnit(float temperature_K, NTC_TemperatureUnit_t unit){
    switch(unit){
        case NTC_MILLI_KELVIN:
            return (int32_t) (temperature_K * 1000.0);
            
        case NTC_MILLI_DEG_CELSIUS:
            return (int32_t) ((temperature_K - 273.15) * 1000.0);
            
        //TODO: implement conversion to Fahrenheit
            
        default: 
            return 0;
    }
}

static float NTC_unitToKelvin(int32_t temperature, NTC_TemperatureUnit_t unit){
    //convert from a temperature value with specified unit to Kelvin
    switch(unit){
        case NTC_MILLI_KELVIN:
            return (float) temperature / 1000.0;
            
        case NTC_MILLI_DEG_CELSIUS:
            return ((float) temperature / 1000.0) + 273.15;
            
        //TODO: implement conversion from Fahrenheit
            
        default: 
            return 0;
    }
}















/*
 * Config file tool. Looks through a given file and finds the first line that contains a definition of it
 * 
 * Line format: 
 *      //a very nice comment, as long as it contains "//" as the first two non space characters
 *      aVeryNiceKeyName = a very nice value with spaces that ends with a newline char //a comment after this is also possible :)
 * 
 *      WARNING: maximum line length given by CONFIG_MAX_LINE_SIZE, entries that are longer will be ignored
 *      WARNING: maximum line count given by CONFIG_MAX_LINE_COUNT, no more lines will be read
 *      NOTE: value will be trimmed of leading and lagging spaces
 */

#if !__is_compiling || __has_include("ff.h")
char * CONFIG_getKey(FIL * file, char * keyToFind){
    //check if the file is value
    if(file < 0xff){
        //no its not, either NULL or an error code. Just return
        return NULL;
    }
    
    //file is valid, seek to the beginning
    FRESULT seekRes = f_lseek(file, 0);
    if(seekRes != FR_OK) TERM_printDebug(TERM_handle, "seek failed (%d)", seekRes);
    
    uint32_t rowsRead = 0;
    char * rowBuffer = pvPortMalloc(CONFIG_MAX_LINE_SIZE * sizeof(char));
    if(rowBuffer == NULL) return NULL;  //no space available for row buffer
    
    //go through the file
    while(rowsRead < CONFIG_MAX_LINE_COUNT){
        //read a line
        uint32_t bytesRead = f_gets(rowBuffer, CONFIG_MAX_LINE_SIZE * sizeof(char), file);
        
        //did it work?
        if(bytesRead == 0){
            //no, some error occurred before we found the key, just return null after cleaning up
            //TERM_printDebug(TERM_handle, "didn't manage to read even one byte :( ");
            
            //free the row buffer
            vPortFree(rowBuffer);        
            return NULL;
            
        }else if(bytesRead == CONFIG_MAX_LINE_SIZE * sizeof(char)){
            //we came across a line too big to fit our buffer, try to find the end of this line so we can skip it properly
            while((bytesRead = f_gets(rowBuffer, CONFIG_MAX_LINE_SIZE * sizeof(char), file))){ 
                if(bytesRead < CONFIG_MAX_LINE_SIZE * sizeof(char)) break;
                //TODO check if we somehow overran the buffer
            }
            //TERM_printDebug(TERM_handle, "line[%d] was too large an skipped ", rowsRead);
            
            //continue scanning the next line
            rowsRead++;
            continue;
            
        }else if(bytesRead > CONFIG_MAX_LINE_SIZE * sizeof(char)){
            //we overran the buffer!! code is now in an instable state, reboot is unfortunately the only option now
            //TODO
            while(1);   //wait for reboot
        }else if(bytesRead <= 2){
            //empty line :) just skip it
            continue;
        }
        
        //TERM_printDebug(TERM_handle, "got line[%d] with length %d: \"%s\" ", rowsRead, bytesRead, rowBuffer);
        
        enum {key_trimLeadingSpaces, key_findEnd, equalSign_find, value_trimLeadingSpaces, value_findEnd} state = key_trimLeadingSpaces;
        
        //check the current line for two patterns: "//" and "="
        uint32_t consequitiveslashCount = 0;
        
        char * key = NULL;
        uint32_t keyEndFound = 0;
        
        uint32_t equalSignFound = 0;
        
        char * value = NULL;
        char * currentEndOfValue = NULL;
        char * currentEndOfValueBeforeCommentSequence = NULL;
        uint32_t wasWaitingForValueBeforeCommentSequence = 0;
        
        for(uint32_t currChar = 0; currChar < bytesRead-1; currChar++){
            //check what letter we are scanning
            
            //TERM_printDebug(TERM_handle, "\tscanning letter '%c' (%02x) ", rowBuffer[currChar], rowBuffer[currChar]);
            
            //are we at the end of the string?
            if(rowBuffer[currChar] == 0) break; //yes => exit loop
            
            //check if we have run across a comment sequence ("//")
            if(rowBuffer[currChar] == '/'){
                //yes, now check if the previous char was a / as well
                if(consequitiveslashCount == 1){
                    //yes, we just scanned across a comment start sequence. 
                    
                    //TERM_printDebug(TERM_handle, "\t\tend of command sequence ");
                    
                    if(currentEndOfValueBeforeCommentSequence != NULL){ 
                        currentEndOfValue = currentEndOfValueBeforeCommentSequence;
                        //TERM_printDebug(TERM_handle, "\t\twe were waiting for a non space char before this comment sequence ");
                    }
                    
                    if(wasWaitingForValueBeforeCommentSequence){
                        value = NULL;
                        state = value_trimLeadingSpaces;
                        //TERM_printDebug(TERM_handle, "\t\twe were waiting for a non space char after equal sign before this comment sequence ");
                    }
                    
                    //cut the string off before the comment sequence and stop scanning
                    rowBuffer[currChar-1] = 0;
                    break;
                }else{
                    //TERM_printDebug(TERM_handle, "\t\tstart of command sequence ");
                    consequitiveslashCount++;   //keep track of how many slashes in a row we found
                    
                    if(currentEndOfValue != NULL) currentEndOfValueBeforeCommentSequence = currentEndOfValue;
                    if(state == value_trimLeadingSpaces) wasWaitingForValueBeforeCommentSequence = 1;
                }
            }else{
                consequitiveslashCount = 0;
                wasWaitingForValueBeforeCommentSequence = 0;
                currentEndOfValueBeforeCommentSequence = NULL;
            }
            
            //run detection statemachine
            switch(state){
                case key_trimLeadingSpaces:
                    //we haven't run across any non space characters so far, check if this is one
                    if(rowBuffer[currChar] != ' ' && !isAsciiSpecialCharacter(rowBuffer[currChar])){
                        //yay start of potential key found :)
                        
                        //now make sure this isn't the key-to-value seperator already
                        if(rowBuffer[currChar] == '='){
                            //hmmm it is, that means no key was found and we may as well stop analysing this line
                            //TERM_printDebug(TERM_handle, "\t\tfirst non space char after start of string, but its '=' so skipping this line ");
                            currChar = bytesRead; //set current pointer to something outside the bounds of the array
                            break;
                        }
                        
                        //TERM_printDebug(TERM_handle, "\t\tfirst non space char after start of string ");
                        
                        //yep a valid first char of the key, note down the position
                        key = &rowBuffer[currChar];
                        
                        //switch to next state: finding the end of the key (first non letter char or the equal sign)
                        state = key_findEnd;
                    }
                    break;
                    
                case key_findEnd:
                    //is this a space that would signal the end of the key?
                    if(isAsciiSpecialCharacter(rowBuffer[currChar]) || rowBuffer[currChar] == ' '){
                        //TERM_printDebug(TERM_handle, "\t\tfirst non letter after start of key ");
                        //yes :) add a string terminator here
                        rowBuffer[currChar] = 0;
                        
                        //switch to next state: finding the equal sign
                        state = equalSign_find;
                    }else if(rowBuffer[currChar] == '='){
                        //TERM_printDebug(TERM_handle, "\t\tfirst non letter after start of key and an equal sign too ");
                        //oh wow already found the equal sign, that also marks the end of the key string
                        rowBuffer[currChar] = 0;
                        
                        //skip next state and go straight to finding the start of the value string
                        state = value_trimLeadingSpaces;
                    }
                    break;
                    
                case equalSign_find:
                    if(rowBuffer[currChar] == '='){
                        //TERM_printDebug(TERM_handle, "\t\tequal sign too ");
                        //found the equal sign
                        
                        //switch to next state: finding the start of the value string
                        state = value_trimLeadingSpaces;
                    }
                    break;
                    
                case value_trimLeadingSpaces:
                    //we haven't run across any non space characters after the equal sign so far, check if this is one
                    if(rowBuffer[currChar] != ' ' && !isAsciiSpecialCharacter(rowBuffer[currChar])){
                        //TERM_printDebug(TERM_handle, "\t\tfirst non space char after equal sign, start of value ");
                        //yay start of potential value string found :)
                        
                        //note down the position
                        value = &rowBuffer[currChar];
                        
                        //switch to next state: finding the end of the key (first non letter char or the equal sign)
                        state = value_findEnd;
                    }
                    break;
                    
                case value_findEnd:
                    //is this a space that we might need to trim?
                    
                    if(isAsciiSpecialCharacter(rowBuffer[currChar]) || rowBuffer[currChar] == ' '){
                        //yes, check if we have already got a potential candidate for the end of the value string
                        if(currentEndOfValue == NULL){
                            //TERM_printDebug(TERM_handle, "\t\tpotential first space after value string ");
                            
                            //no, no candidate exists. That means we are a suspect (pretty sus if you ask me)
                            currentEndOfValue = &rowBuffer[currChar];
                        }else; //there is already another candidate. Aka no non-space char was scanned since entering the current space region 
                        
                    }else{
                        //ah we got a non space char, reset the potential end of the buffer
                        if(currentEndOfValue != NULL){
                            currentEndOfValue = NULL;
                            //TERM_printDebug(TERM_handle, "\t\tnope*d ");
                        }
                    }
                    break;
            }
        }
        
        //row scanning statemachine has exited, see how far it got
        if(state == value_findEnd){
            //is there a trailing space char we need to trim?
            if(currentEndOfValue != NULL) *currentEndOfValue = 0; //yes, make a terminator out of it (t800 should do fine for now, might upgrade to 1000 in the future/past)
                
            //got all the way through to the end, so a potential key and value pair exists in the buffer now
            //TERM_printDebug(TERM_handle, "\tstate machine got all the way through :D ");
            //TERM_printDebug(TERM_handle, "\t\t key=\"%s\" value=\"%s\" ", key, value);
            
            //now check if the key matches the one we are looking for
            if(strcmp(keyToFind, key) == 0){
                //yes we found it!! :D now transfer the value string into its own buffer and return it
                
                //create and copy to new buffer
                char * ret = pvPortMalloc(strlen(value)+1);
                strcpy(ret, value);
                vPortFree(rowBuffer);
                
                TERM_printDebug(TERM_handle, "key matches! final value return=\"%s\" ", ret);
                return ret;
            }
        }
                
        //TERM_printDebug(TERM_handle, "row wasn't what we hoped :( \r\n\n\n\n");
        
        rowsRead++;
    }
    
    //TERM_printDebug(TERM_handle, "scanned %d rows and didn't find the key :(\r\n", rowsRead);
        
    //free the row buffer
    vPortFree(rowBuffer);   
    return NULL;
}
#endif

uint32_t isAsciiSpecialCharacter(char c){
    if(c < 32) return 1;
    if(c > 127) return 1;
    return 0;
}

uint32_t isAsciiNumber(char c){
    if(c >= '0' && c <= '9') return 1;
    return 0;
}






/*
 * Fixed point atoi function that supports multipliers like k,m,M,c,d etc.
 * 
 * input string must start with a number, otherwise 0 is returned. It will scan until either a null terminator is reached or strlen chars have been read
 * 
 * ignoreUnit will make the function not return in an error state if the string contains a unit at the end (like for example "1sec")
 * 
 *      NOTE: if the returned value would contain an integer overflow due to an excessively large input number a 0 is returned instead
 *      NOTE: any numbers after the last digit in the fixed point return will be ignored
 */

static int32_t ctoi(char c){
    if(!isAsciiNumber(c)) return 0;
    return c - '0';  //if we are sure that the char is a number, then we can convert it by subtracting the offset of ascii '0' from it
}

static int32_t getExponent(char c){
    if(c == 'f')    return -15;
    if(c == 'p')    return -12;
    if(c == 'n')    return -9;
    if(c == 'u' || c == 230) return -6;    //ascii 230 is the greek micro
    if(c == 'm')    return -3;
    if(c == 'c')    return -2;
    if(c == 'd')    return -1;
    if(c == '.')    return 0;
    if(c == 'h')    return 2;
    if(c == 'k')    return 3;
    if(c == 'M')    return 6;
    if(c == 'G')    return 9;
    if(c == 'T')    return 12;
    if(c == 'P')    return 15;
    
    if(isAsciiNumber(c)) return 0xffffffff;
    return 0x7fffffff;
}

//precalculated exponents that are reasonable for a int32_t number
static int32_t exponentTable[10] = {
1, //10 ^ 0
10, //10 ^ 1
100, //10 ^ 2
1000, //10 ^ 3
10000, //10 ^ 4
100000, //10 ^ 5
1000000, //10 ^ 6
10000000, //10 ^ 7
100000000, //10 ^ 8
1000000000, //10 ^ 9
};

int32_t atoiFP(char * a, uint32_t strlen, int32_t baseExponent, uint32_t ignoreUnit){
    int32_t ret = 0;
    
    enum {waitStart, waitMultiplier, waitEnd, waitE_Exponent} state = waitStart;
    
    int32_t currentNumberPositionExponent = 0;
    
    uint32_t isNumberNegative = 0;
    uint32_t isExponentADot = 0;
    uint32_t breakBit = 0;
    
    char* startOfNumber = NULL;
    char* startOfNumberExponent = NULL;
    
    uint32_t checkForLength = strlen > 0;
    
    //iterate through the string
    while(*a != 0 && !breakBit){
        if(checkForLength) if(!strlen--) break;
        //TERM_printDebug(TERM_handle, "\r\n\nscanning letter '%c' (%02x) ", *a, *a);
        
        int32_t currentExponent = getExponent(*a);
        
        //which state are we in?
        switch(state){
            case waitStart:
                //skip any leading spaces
                if(*a != ' '){
                    //found a non space char, is it a number?
                    if(isAsciiNumber(*a)){
                        //TERM_printDebug(TERM_handle, "\t found start of number ");
                        //yes, start of number to convert
                        startOfNumber = a;
                        state = waitMultiplier;
                        
                    }else if(*a == '.'){
                        //we got a dot as the first char, skip waiting for a first exponent as we already have one here
                        startOfNumber = a;
                        currentNumberPositionExponent = 0;
                        state = waitEnd;
                        isExponentADot = 1;
                        
                    }else if(*a == '-'){
                        //no, but the current character is a minus sign, invert the sign of the return number
                        isNumberNegative = !isNumberNegative;
                        //TERM_printDebug(TERM_handle, "\tfound sign (number is now %d) ", isNumberNegative ? "negative" : "positive");
                        
                    }else if(*a != '+'){    //a number string with a leading plus is also valid, but we don't need to do anything with it
                        //TERM_printDebug(TERM_handle, "\tfound sign (number is now %d) ", isNumberNegative ? "negative" : "positive");
                        //no, apparently the string doesn't start with a number => invalid :(
                        return 0;
                    }
                }
                break;
            case waitMultiplier:
                if(currentExponent == 0x7fffffff){
                    //current character is neither a valid exponent designator or a number. we likely reached the end of the string we need to parse
                    
                    //check if by any chance this is an 'e' and there is also at least one more char in this string we can scan
                    if((*a == 'e' || *a == 'E') && (!checkForLength || strlen) && *(a+1) != 0){
                        //yes it is
                        if(*(a+1) == '-' || isAsciiNumber(*(a+1))){
                            //TERM_printDebug(TERM_handle, "\tfound start of exponent string ");
                            //and it is also a number, yay. That means we must now scan for the end of the added exponent
                            startOfNumberExponent = a + 1;
                            state = waitE_Exponent;
                        }else{
                            //TERM_printDebug(TERM_handle, "\tpotential start of exponent string was invalid ");
                            //break out of the loop and continue on to the conversion
                            breakBit = 1;
                            break;
                        }
                    }else{
                        //TERM_printDebug(TERM_handle, "\tend of number ");
                        //break out of the loop and continue on to the conversion
                        breakBit = 1;
                        break;
                    }
                    
                }else if(currentExponent != 0xffffffff){
                        
                    //we found an exponent of some description. remember its value as we continue through the string
                    currentNumberPositionExponent = currentExponent;
                    state = waitEnd;
                    
                    //now check if the exponent is a dot, if so we can accept another real exponent at the end of the string
                    if(currentExponent == 0) isExponentADot = 1;
                    
                    //TERM_printDebug(TERM_handle, "\tfound a first exponent (%d) currNumExp=%d ", currentExponent, currentNumberPositionExponent);
                }else{
                    //TERM_printDebug(TERM_handle, "\tthat's just a digit :) ");
                }
                break;
            case waitEnd:
                if(!isAsciiNumber(*a)){
                    //we've reached another non number
                    
                    //check if the first non number found was a dot and we currently have a valid exponent
                    if(isExponentADot){
                        
                        //yes it was. Now check if we have reached the start of a numerical exponent string
                        if((*a == 'e' || *a == 'E') && strlen){
                            //yes it is
                            if(*(a+1) == '-' || isAsciiNumber(*(a+1))){
                                //TERM_printDebug(TERM_handle, "\tfound exponent string start ");
                                //and it is also a number, yay. That means we must now scan for the end of the added exponent
                                startOfNumberExponent = a + 1;
                                state = waitE_Exponent;
                            }else{
                                //TERM_printDebug(TERM_handle, "\tfound non exponent string start ");
                                
                                //break out of the loop and continue on to the conversion
                                breakBit = 1;
                                break;
                            }
                        }else if(currentExponent != 0x7fffffff){
                            //TERM_printDebug(TERM_handle, "\tfound exponent after string %d ", currentExponent);
                            
                            currentNumberPositionExponent += currentExponent;
                    
                            //break out of the loop and continue on to the conversion
                            breakBit = 1;
                            break;
                        }else{
                            //TERM_printDebug(TERM_handle, "\tfound end of number string without exponent after dot");
                            
                            //a dot was the first exponent and we have now reached the first non digit after it, but it isn't a valid exponent (so probably a unit)
                            //check if we are allowed to ignore this and convert anyway or must return with an error
                            if(!ignoreUnit) return 0;

                            //otherwise break out of the loop and continue on to the conversion
                            breakBit = 1;
                            break;
                        }
                    }else{
                        //TERM_printDebug(TERM_handle, "\tfound end of number string with exponent ");
                        //no we already got a valid exponent or aren't even scanning one now. We must be at the start of the unit
                        //check if we are allowed to ignore this and convert anyway or must return with an error
                        if(!ignoreUnit) return 0;
                    
                        //otherwise break out of the loop and continue on to the conversion
                        breakBit = 1;
                        break;
                    }
                }else{
                    //we are scanning a number, make sure to decrease the current numbers exponent as we are after the point
                    currentNumberPositionExponent--;
                    //TERM_printDebug(TERM_handle, "\tcurrNumExp=%d ", currentNumberPositionExponent);
                }
                break;
                
            case waitE_Exponent:
                //we have passed a "e" in the string and must now find the end of the number after it
                
                //check if we are still scanning a number or the sign (but sign only if it is the first char in the exponent string)
                if(!isAsciiNumber(*a) && !((*a == '-' || *a == '+') && a == startOfNumberExponent)){
                    //TERM_printDebug(TERM_handle, "\tfound end of exponent string ");
                    //found first non sign and non number char, that means we are now at the end of the exponent number string
                    
                    //break the loop and let the conversion routinge handle the exponent string
                    breakBit = 1;
                    break;
                }else{
                    //TERM_printDebug(TERM_handle, "\tthat's just a digit :) ");
                }
                break;
        }
        
        if(!breakBit) a++;
    }
    
    //are we still waiting for the end on an exponent string?
    if(state == waitE_Exponent){
        //yes, we are now at the end of it
        //TERM_printDebug(TERM_handle, "string ended with end of exponent string ");
        //found first non sign and non number char, that means we are now at the end of the exponent number string

        //copy that number into a new string for conversion
        uint32_t i = 0;
        char exponentBuffer[11];
        char * exponentStringPosition = startOfNumberExponent;
        for(; i < 10 && exponentStringPosition != a; i++){
            exponentBuffer[i] = *(exponentStringPosition++);
        }
        //terminate the string
        exponentBuffer[i] = 0;

        //did we actually get an entire string or was the buffer too small?
        if(exponentStringPosition == a){
            //yes that worked, now convert the number to an int with the normal atoi function
            int32_t newExponent = atoi(exponentBuffer);
            currentNumberPositionExponent += newExponent;

            //TERM_printDebug(TERM_handle, "\t\tcopied string: \"%s\"=", exponentBuffer, newExponent);

            //now roll back the current char pointer to the 'e' we scanned before entering this state
            //this is to allow for the conversion to take place as if no numerical exponent existed after the number
            a = startOfNumberExponent-1;
        }else{
            //TERM_printDebug(TERM_handle, "\t\texponent buffer was too small :(");
            //no didn't work unfortunately, that means the number read was invalid and we must return 0
            return 0;
        }
    }
    
    //check if we successfully exited the loop
    if(startOfNumber == NULL){
        //no :( no number was found
        return 0;
    }
    
    //shift the fixed point exponent by the base exponent
    currentNumberPositionExponent += baseExponent;
                        
    //TERM_printDebug(TERM_handle, "Now converting the number with the current exponent of %d (from base %d)", currentNumberPositionExponent, baseExponent);
    
    //start of number found, we are now one byte after the end of the number string, step backwards through the string and do the converstion
    while(*--a != 0 && a >= startOfNumber){
        //TERM_printDebug(TERM_handle, "\r\n\nscanning letter '%c' (%02x) ", *a, *a);
        
        //is the character we are at right now a number?
        if(isAsciiNumber(*a)){
            //yes, would it actually be put into the return value? it wouldn't if the exponent was too low or too big
            if(currentNumberPositionExponent >= 0){
                if(currentNumberPositionExponent < 10){
                    //yes, convert it
                    int32_t currentNumber = ctoi(*a);
                    
                    //TERM_printDebug(TERM_handle, "\tconverted digit to int = %d", currentNumber);
                    
                    //and finally do a check for a potential overflow if a digit larger than 2 gets multiplied with the largest exponent
                    if(currentNumberPositionExponent == 9 && currentNumber > 2){
                        //TERM_printDebug(TERM_handle, "\t=> would overflow :( ");
                        //yes overflow would occur :(
                        return 0;
                    }
                    
                    //no overflow would occur, scale the number by the current exponent
                    currentNumber *= exponentTable[currentNumberPositionExponent];
                    
                    //TERM_printDebug(TERM_handle, "\tscaled with exponent %d (10^%d) = %d ", exponentTable[currentNumberPositionExponent], currentNumberPositionExponent, currentNumber);
                    
                    ret += currentNumber;
                    
                    //TERM_printDebug(TERM_handle, "\tret is now %d ", ret);
                }else{
                    //TERM_printDebug(TERM_handle, "\tignoring number due to exponent > 10 (%d)", currentNumberPositionExponent);
                    //a number with too big of an exponent would be converted, which would cause an overflow. return 0 so the user notices something is wrong
                    return 0;
                }
            }else{
                //TERM_printDebug(TERM_handle, "\tignoring number due to exponent < 0 (%d)", currentNumberPositionExponent);
            }
                    
            //increase position exponent by one for the digit to the left of us
            currentNumberPositionExponent++;
        }else{
            //TERM_printDebug(TERM_handle, "\tnot a number, skipping character ");
            //current char is not a number, this should only occur once [perhaps TODO check if that has worked] (as we scan over the point or the exponent). Just skip it
        }
    }
                    
    //TERM_printDebug(TERM_handle, "conversion done! result without sign = %d ", ret);
    
    //and finally return the number taking the sign bit into account
    return isNumberNegative ? -ret : ret;
}









/*
 * LUT based fast sine function
 * 
 * parameter is a variable between 0-512 with 512 = 2pi. Return value is 0-500000
 *     
 */


//256 value sin lookup table x := [0, pi)
//TODO reduce this by half, symmetry's a bitch ;)
static const int32_t sineTable[] = {0, 6135, 12270, 18403, 24533, 30660, 36782, 42898, 49008, 55111, 61205, 67290, 73365, 79429, 85480, 91519, 97545, 103555, 109550, 115529, 121490, 127432, 133356, 139259, 145142, 151002, 156840, 162655, 168444, 174209, 179947, 185658, 191341, 196996, 202620, 208214, 213777, 219308, 224805, 230269, 235698, 241091, 246449, 251769, 257051, 262294, 267498, 272662, 277785, 282865, 287904, 292898, 297849, 302755, 307615, 312429, 317196, 321915, 326586, 331207, 335779, 340300, 344770, 349188, 353553, 357865, 362123, 366327, 370475, 374568, 378604, 382583, 386505, 390368, 394173, 397918, 401603, 405228, 408792, 412294, 415734, 419112, 422426, 425677, 428864, 431986, 435043, 438035, 440960, 443819, 446612, 449337, 451994, 454583, 457104, 459556, 461939, 464253, 466496, 468669, 470772, 472803, 474764, 476653, 478470, 480215, 481888, 483488, 485015, 486469, 487851, 489158, 490392, 491552, 492638, 493650, 494588, 495451, 496239, 496953, 497592, 498156, 498645, 499059, 499397, 499661, 499849, 499962, 500000, 499962, 499849, 499661, 499397, 499059, 498645, 498156, 497592, 496953, 496239, 495451, 494588, 493650, 492638, 491552, 490392, 489158, 487851, 486469, 485015, 483488, 481888, 480215, 478470, 476653, 474764, 472803, 470772, 468669, 466496, 464253, 461939, 459556, 457104, 454583, 451994, 449337, 446612, 443819, 440960, 438035, 435043, 431986, 428864, 425677, 422426, 419112, 415734, 412294, 408792, 405228, 401603, 397918, 394173, 390368, 386505, 382583, 378604, 374568, 370475, 366327, 362123, 357865, 353553, 349188, 344770, 340300, 335779, 331207, 326586, 321915, 317196, 312429, 307615, 302755, 297849, 292898, 287904, 282865, 277785, 272662, 267498, 262294, 257051, 251769, 246449, 241091, 235698, 230269, 224805, 219308, 213777, 208214, 202620, 196996, 191341, 185658, 179947, 174209, 168444, 162655, 156840, 151002, 145142, 139259, 133356, 127432, 121490, 115529, 109550, 103555, 97545, 91519, 85480, 79429, 73365, 67290, 61205, 55111, 49008, 42898, 36782, 30660, 24533, 18403, 12270, 6135};

int32_t qSin(int32_t x){
    int32_t xt = abs(x) % 0xff;
    
    if(xt == 0) return 0;
    return (x > 0) ? sineTable[xt] : -sineTable[xt];
}






const uint32_t util_bitMasks[] = {0x00000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000};


/*
 * Bit and byte manipulation 
 */

//gets a 16bit word from an array of bytes
//WARNING: user must make sure that offset+1 bytes are available in the data buffer!
uint16_t get16BitWord(uint8_t * data, uint32_t offset, Endiannes_t endian){
    if(endian == LITTLE_ENDIAN) return (data[offset] << 8) | data[offset + 1];
    if(endian == BIG_ENDIAN) return (data[offset+1] << 8) | data[offset];
}

//sets a 16bit word in an array of bytes
//WARNING: user must make sure that the data buffer is large enough!
uint16_t set16BitWord(uint8_t * data, uint32_t offset, Endiannes_t endian, uint32_t value){
    if(endian == LITTLE_ENDIAN){ 
       data[offset] = (value >> 8) & 0xff;
       data[offset + 1] = value & 0xff;
    }
    
    if(endian == BIG_ENDIAN){ 
       data[offset + 1] = (value >> 8) & 0xff;
       data[offset] = value & 0xff;
    }
}


//get a bit out of a number of bits in a byte array
//WARNING: user must make sure that enough bits are in the data buffer!
uint32_t getBit(uint8_t * data, uint32_t bitNumber){
    //get byte offset
    uint32_t byteOffset = bitNumber / 8;
    uint32_t bitOffset = bitNumber - (byteOffset * 8);
    
    //return 1 if the required bit in the target byte is set
    return (data[byteOffset] & util_bitMasks[bitOffset]) > 1;
}

//set a bit in a byte array
//WARNING: user must make sure that the data buffer is large enough to even have the bit in question!
void setBit(uint8_t * data, uint32_t bitNumber, uint32_t value){
    //get byte offset
    uint32_t byteOffset = bitNumber / 8;
    uint32_t bitOffset = bitNumber - (byteOffset * 8);
    
    //apply the operation to the byte in question
    if(value) data[byteOffset] |= util_bitMasks[bitOffset]; else data[byteOffset] &= ~util_bitMasks[bitOffset];
}

int32_t ceil_div(int32_t x, int32_t y) {
    return (x + y - 1) / y;
}