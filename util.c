#include <stdint.h>
#include <string.h>

#include "util.h"
#include "ff.h"
#include "TTerm.h"
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
 *      {xValue,yValue,dX/dY}
 *      - x/yValue: as in preComputedDerivative=0
 *      - dX/dY: pre-computed rate of change in between the point and the next one in the list
 *  
 *      NOTE: the list must be sorted by x values in ascending order (x[0] < x[1] < x[2]...)
 *      NOTE: if you need this to be fast, make sure to pre-calculate the derivatives as this saves a division for every conversion
 */
int32_t PWL_getY(int32_t x, int32_t * pwl, uint32_t listSizeRows, uint32_t preComputedDerivative){
    if(listSizeRows < 2){
        //can't approximate any function if all we have is a single point. We need at least two
        //configASSERT(0);
        return 0;
    }
    
    int32_t * currentRow = NULL;
    int32_t * lastRow = NULL;
    
    //find the neighbouring points
    int32_t entry = 0;
    for(entry = 0; entry < listSizeRows; entry++){
        //calculate data pointer to current row. If preComputedDiff is on then there is a third value in the row (dX/dY between the current and next entry) otherwise only 2
        currentRow = pwl + entry * (preComputedDerivative ? 2 : 3);
        
        //check if the entry we are looking at has a x value larger than the one we are looking for
        if(currentRow[0] > x){
            //entry found :)
            
            //check if we are at the first entry in the list
            if(lastRow == NULL){ 
                //yes we are, that means that the x value we are looking for lies outside the boundary of the pwl. Try to interpolate with the first two points
                
                //is there at least one point ahead of us?
                if(entry + 1 >= listSizeRows){
                    //hmm no, this shouldn't technicall be possible as we already check for a list size of at least 2. just return, we can't do any calculations now anyway
                    
                    //this btw means that we are left of the first point in the graph, but no other point exists after the first one (aka we only have one in total)
                    return 0;
                }
                
                //change the pointers to point to the correct points
                lastRow     = currentRow;
                currentRow  = pwl + (entry + 1) * (preComputedDerivative ? 2 : 3);
                
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
        //this means loop exited due to entry == listSizeRows. Move lastRow back one entry
        
        //can we actually do that?
        if(entry - 2 < 0){
            //hmm no, this again means there is only one point in the list... can't do anything
            return 0;
        }
        
        lastRow = pwl + (entry-2) * (preComputedDerivative ? 2 : 3);
    }
    
    //now interpolate between the two points pointed to by lastRow (start point, left of X) and currentRow (end point, right of x))
            
    //first calculate the x offset from the start point
    int32_t localX = x - lastRow[0];
    int32_t localY = lastRow[1];        

    //calculate derivative of pwl between start and end points
    int32_t dYdX = 0;
    if(preComputedDerivative){
        //nothing to compute, just read the value from the start point dataset
        dYdX = lastRow[2];
    }else{
        //dy = currentRow[1] - lastRow[1];
        //dx = currentRow[0] - lastRow[0];
        dYdX = (currentRow[1] - lastRow[1]) / (currentRow[0] - lastRow[0]);
    }

    //calculate linear function as 
    //   y = m     *   x       + b         with m=dYdX, x=localX, b=localY
    return   dYdX  *   localX  + localY;
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
                
                //TERM_printDebug(TERM_handle, "key matches! final value return=\"%s\" ", ret);
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

uint32_t isAsciiSpecialCharacter(char c){
    if(c < 32) return 1;
    if(c > 127) return 1;
}