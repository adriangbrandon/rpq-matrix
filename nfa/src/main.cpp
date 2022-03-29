
#define NULL 0
typedef int bool;


bool removeSpaces(char *inStr, char* outStr)
{
    int inPos = 0;
    int outputPos = 0;

    if(inStr == NULL || outStr == NULL)
        return 0;

    //por cada caracter en la posicion inPos
    while (inStr[inPos] > 0)
    {
        // si no es un espacio lo muevo a outPos
        if (inStr[inPos] != ' ')
        {
            if (outputPos < inPos)
                outStr[outputPos] = inStr[inPos];
            ++outputPos;
        }

        ++inPos;
    }

    outStr[outputPos] = inStr[inPos];

    return 1;
};

// inPos: 0
// outPos: 0
//  inStr[inPos] = ' '
//    inPos = 1
// inPos: 1
//  inStr[inPos] = 'h'
//  str: 'hh i \0'
//  outPos: 1

// inPos: 2
// outPos: 1
//  inStr[inPos] = ' '
//  str: 'hh i \0'

// inPos: 3
// outPos: 1
//  inStr[inPos] = 'i'
//         | | 
//  str: 'hh i \0'
//  -> str: 'hi i \0'
//  -> outputPos: 2


// inPos: 4
// outPos: 2
//  inStr[inPos] = ' '
//  str: 'hi i \0'

// inPos: 5
// outPos: 2
//  inStr[inPos] = ' '
//  str: 'hi\0i \0'
