/*Copyright (C) 2018 Wayne Mogg All rights reserved.

This file may be used either under the terms of:

1. The GNU General Public License version 3 or higher, as published by
the Free Software Foundation, or

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef windowedops_h
#define windowedops_h

/*+
________________________________________________________________________

 Author:        Wayne Mogg
 Date:          December 2018
 ________________________________________________________________________

-*/ 
#include "arrayndimpl.h"

namespace windowedOps{
    
void sum( const Array1DImpl<float>& input, const int winSize, Array1DImpl<float>& output )
{
    int sz = input.info().getSize(0);
    output.setSize(sz);
    
    if (sz<winSize) {
        double tmp = 0.0;
        for (int idx=0; idx<sz; idx++)
            tmp += input.get(idx);
        output.setAll(tmp);
    } else {
        double tmp = 0.0;
        int halfWinSize = (winSize-1)/2;
        for (int idx=0; idx<=halfWinSize; idx++)
            tmp += input.get(idx);
        output.set(0, tmp);
        
        for (int idx=1; idx<sz; idx++) {
            int isub = idx - halfWinSize -1;
            double sub = isub<0 ? 0.0 : input.get(isub);
            int iadd = idx + halfWinSize;
            double add = iadd>=sz ? 0.0 : input.get(iadd);
            tmp += add-sub;
            output.set( idx, tmp );
        }
    }
} 

void min( const Array1DImpl<float>& input, const int winSize, Array1DImpl<float>& output )
{
    int sz = input.info().getSize(0);
    output.setSize(sz);

    if (sz<winSize) {
        float tmp = input.get(0);
        for (int idx=0; idx<sz; idx++)
            tmp = tmp<=input.get(idx) ? tmp: input.get(idx);
        output.setAll(tmp);
    } else {
        float tmp = input.get(0);
        int halfWinSize = (winSize-1)/2;
        for (int idx=0; idx<=halfWinSize; idx++)
            tmp = tmp<=input.get(idx) ? tmp: input.get(idx);
        output.set(0, tmp);
        
        for (int idx=1; idx<sz; idx++) {
            int isub = idx - halfWinSize -1;
            isub = isub<0? 0: isub;
            int iadd = idx + halfWinSize;
            iadd = iadd<sz ? iadd: sz-1;  
            if (tmp == input.get(isub)) {
                tmp = input.get(iadd);
                for (int i=iadd; i>isub; i--)
                    tmp = tmp<=input.get(i) ? tmp: input.get(i);
            }
            output.set( idx, tmp);
        }
    }
} 

void max( const Array1DImpl<float>& input, const int winSize, Array1DImpl<float>& output )
{
    int sz = input.info().getSize(0);
    output.setSize(sz);

    if (sz<winSize) {
        float tmp = input.get(0);
        for (int idx=0; idx<sz; idx++)
            tmp = tmp>=input.get(idx) ? tmp: input.get(idx);
        output.setAll(tmp);
    } else {
        float tmp = input.get(0);
        int halfWinSize = (winSize-1)/2;
        for (int idx=0; idx<=halfWinSize; idx++)
            tmp = tmp>=input.get(idx) ? tmp: input.get(idx);
        output.set(0, tmp);
        
        for (int idx=1; idx<sz; idx++) {
            int isub = idx - halfWinSize -1;
            isub = isub<0? 0: isub;
            int iadd = idx + halfWinSize;
            iadd = iadd<sz ? iadd: sz-1;  
            if (tmp == input.get(isub)) {
                tmp = input.get(iadd);
                for (int i=iadd; i>isub; i--)
                    tmp = tmp>=input.get(i) ? tmp: input.get(i);
            }
            output.set( idx, tmp);
        }
    }
} 

void minIdx( const Array1DImpl<float>& input, const int winSize, Array1DImpl<int>& output )
{
    int sz = input.info().getSize(0);
    output.setSize(sz);

    if (sz<winSize) {
        int tmp = 0;
        for (int idx=0; idx<sz; idx++)
            tmp = input.get(tmp)<=input.get(idx) ? tmp: idx;
        output.setAll(tmp);
    } else {
        int tmp = 0;
        int halfWinSize = (winSize-1)/2;
        for (int idx=0; idx<=halfWinSize; idx++)
            tmp = input.get(tmp)<=input.get(idx) ? tmp: idx;
        output.set(0, tmp);
        
        for (int idx=1; idx<sz; idx++) {
            int isub = idx - halfWinSize -1;
            isub = isub<0? 0: isub;
            int iadd = idx + halfWinSize;
            iadd = iadd<sz ? iadd: sz-1;  
            if (tmp == isub) {
                tmp = iadd;
                for (int i=iadd; i>isub; i--)
                    tmp = input.get(tmp)<=input.get(i) ? tmp: i;
            }
            output.set( idx, tmp);
        }
    }
} 

void maxIdx( const Array1DImpl<float>& input, const int winSize, Array1DImpl<int>& output )
{
    int sz = input.info().getSize(0);
    output.setSize(sz);

    if (sz<winSize) {
        int tmp = input.get(0);
        for (int idx=0; idx<sz; idx++)
            tmp = input.get(tmp)>=input.get(idx) ? tmp: idx;
        output.setAll(tmp);
    } else {
        int tmp = input.get(0);
        int halfWinSize = (winSize-1)/2;
        for (int idx=0; idx<=halfWinSize; idx++)
            tmp = input.get(tmp)>=input.get(idx) ? tmp: idx;
        output.set(0, tmp);
        
        for (int idx=1; idx<sz; idx++) {
            int isub = idx - halfWinSize -1;
            isub = isub<0? 0: isub;
            int iadd = idx + halfWinSize;
            iadd = iadd<sz ? iadd: sz-1;  
            if (tmp == isub) {
                tmp = iadd;
                for (int i=iadd; i>isub; i--)
                    tmp = input.get(tmp)>=input.get(i) ? tmp: i;
            }
            output.set( idx, tmp);
        }
    }
} 

}
#endif
