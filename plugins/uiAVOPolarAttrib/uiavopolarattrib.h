/*Copyright (C) 2018 Wayne Mogg. All rights reserved.

This file may be used either under the terms of:

1. The GNU General Public License version 3 or higher, as published by
the Free Software Foundation, or

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef uiavopolarattrib_h
#define uiavopolarattrib_h

/*+
________________________________________________________________________

 Author:        Wayne Mogg
 Date:          November 2018
 ________________________________________________________________________

-*/

#include "uiavopolarattribmod.h"
#include "uiattrdesced.h"

class uiAttrSel;
class uiGenInput;
class uiStepOutSel;


/*! \brief Compute AVO Polarization, Strength and other related attributes description editor */

class uiAVOPolarAttrib : public uiAttrDescEd
{ mODTextTranslationClass(uiAVOPolarAttrib);
public:

    uiAVOPolarAttrib(uiParent*,bool);

protected:

    uiAttrSel*		inp_interceptfld_;
    uiAttrSel*		inp_gradientfld_;
    uiGenInput*		outputfld_;
    uiGenInput*     gateBGfld_;
    uiStepOutSel*   stepoutBGfld_;
    uiGenInput*     gatefld_;

    bool        setParameters(const Attrib::Desc&);
    bool        setInput(const Attrib::Desc&);
    bool		setOutput(const Attrib::Desc&);
    bool        getParameters(Attrib::Desc&);
    bool        getInput(Attrib::Desc&);
    bool		getOutput(Attrib::Desc&);

                mDeclReqAttribUIFns
};


#endif


