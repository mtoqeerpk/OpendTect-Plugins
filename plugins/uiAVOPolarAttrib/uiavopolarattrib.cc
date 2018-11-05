/*Copyright (C) 2018 Wayne Mogg. All rights reserved.

This file may be used either under the terms of:

1. The GNU General Public License version 3 or higher, as published by
the Free Software Foundation, or

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*/

/*+
________________________________________________________________________

 Author:        Wayne Mogg
 Date:          November 2018
 _______________________________________________________________________

-*/

#include "uiavopolarattrib.h"
#include "avopolarattrib.h"

#include "attribdesc.h"
#include "attribparam.h"
#include "uiattribfactory.h"
#include "uiattrsel.h"
#include "uigeninput.h"

using namespace Attrib;

static const char* outputstr[] =
{
    "Background Angle",
    "Local Angle",
    "Angle Difference",
    "Strength",
    "Polarization Product",
    "Quality",
    "Background Quality",
    0
};


mInitAttribUI(uiAVOPolarAttrib,AVOPolarAttrib,"AVO Polarization Attribute",sKeyFilterGrp())


uiAVOPolarAttrib::uiAVOPolarAttrib( uiParent* p, bool is2d )
    : uiAttrDescEd(p,is2d,"mToDoHelpID")

{
    inp_interceptfld_ = createInpFld( is2d, "Intercept Volume" );
    inp_gradientfld_ = createInpFld( is2d, "Gradient Volume" );
    inp_gradientfld_->attach( alignedBelow, inp_interceptfld_ );

    outputfld_ = new uiGenInput( this, uiStrings::sOutput(),StringListInpSpec(outputstr) );
    outputfld_->attach( alignedBelow, inp_gradientfld_ );

    zmarginBGfld_ = new uiGenInput( this, "Background Analysis Z Window (samples)", IntInpIntervalSpec().setName("Samples after",1).setName("Samples before",0) );
    zmarginBGfld_->attach( alignedBelow, outputfld_ );
    zmarginBGfld_->valuechanging.notify(mCB(this, uiAVOPolarAttrib,doZmarginCheck) );

    stepoutBGfld_ = new uiStepOutSel( this, is2d );
    stepoutBGfld_->attach( rightTo, zmarginBGfld_ );
    stepoutBGfld_->valueChanging.notify(mCB(this, uiAVOPolarAttrib,doStepOutCheck) );
    stepoutBGfld_->setFieldNames( "Inl Stepout", "Crl Stepout" );

    zmarginfld_ = new uiGenInput( this, "Local Analysis Z Window (samples)", IntInpIntervalSpec().setName("Samples after",1).setName("Samples before",0) );
    zmarginfld_->attach( alignedBelow, zmarginBGfld_ );
    zmarginfld_->valuechanging.notify(mCB(this, uiAVOPolarAttrib,doZmarginCheck) );

    
    setHAlignObj( outfld_ );
}




bool uiGradientAttrib::setParameters( const Attrib::Desc& desc )
{
    if ( desc.attribName() != GradientAttrib::attribName() )
	return false;

    mIfGetEnum(GradientAttrib::outputStr(),output, outfld_->setValue(output) )
    mIfGetEnum(GradientAttrib::operatorStr(), opert, operatorfld_->setValue(opert))
    return true;
}


bool uiGradientAttrib::setInput( const Attrib::Desc& desc )
{
    putInp( inpfld_, desc, 0 );
    return true;
}


bool uiGradientAttrib::getParameters( Attrib::Desc& desc )
{
    if ( desc.attribName() != GradientAttrib::attribName() )
	return false;

    mSetEnum( GradientAttrib::outputStr(), outfld_->getIntValue() );
    mSetEnum( GradientAttrib::operatorStr(), operatorfld_->getIntValue() );
    
    return true;
}


bool uiGradientAttrib::getInput( Attrib::Desc& desc )
{
    inpfld_->processInput();
    fillInp( inpfld_, desc, 0 );
    return true;
}



