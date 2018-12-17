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
#include "uistepoutsel.h"
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

    gateBGfld_ = new uiGenInput( this, zDepLabel(tr("Background "), tr("gate")), FloatInpIntervalSpec().setName("Z start",0).setName("Z stop",1));
    gateBGfld_->attach( alignedBelow, outputfld_ );

    stepoutBGfld_ = new uiStepOutSel( this, is2d );
    stepoutBGfld_->attach( rightTo, gateBGfld_ );
    stepoutBGfld_->setFieldNames( is2d ? "Trace Nr Stepout" : "Inl Stepout", "Crl Stepout" );

    gatefld_ = new uiGenInput( this, zDepLabel(tr("Local "), tr("gate")), FloatInpIntervalSpec().setName("Z start",0).setName("Z stop",1));
    gatefld_->attach( alignedBelow, gateBGfld_ );
    
    setHAlignObj( outputfld_ );
}

bool uiAVOPolarAttrib::setParameters( const Attrib::Desc& desc )
{
    if ( desc.attribName() != AVOPolarAttrib::attribName() )
	return false;

    mIfGetFloatInterval(AVOPolarAttrib::gateBGStr(), gateBG, gateBGfld_->setValue(gateBG))
    mIfGetBinID(AVOPolarAttrib::soBGStr(),soBG,stepoutBGfld_->setBinID(soBG))
    mIfGetFloatInterval(AVOPolarAttrib::gateStr(), gate, gatefld_->setValue(gate))
    
    return true;
}


bool uiAVOPolarAttrib::setInput( const Attrib::Desc& desc )
{
    putInp( inp_interceptfld_, desc, 0 );
    putInp( inp_gradientfld_, desc, 1 );
    return true;
}

bool uiAVOPolarAttrib::setOutput( const Desc& desc )
{
    outputfld_->setValue( desc.selectedOutput() );
    return true;
}


bool uiAVOPolarAttrib::getParameters( Attrib::Desc& desc )
{
    if ( desc.attribName() != AVOPolarAttrib::attribName() )
        return false;

    mSetBinID( AVOPolarAttrib::soBGStr(), stepoutBGfld_->getBinID() );
    mSetFloatInterval( AVOPolarAttrib::gateBGStr(), gateBGfld_->getFInterval() );
    mSetFloatInterval( AVOPolarAttrib::gateStr(), gatefld_->getFInterval() );

    return true;
}


bool uiAVOPolarAttrib::getInput( Attrib::Desc& desc )
{
    fillInp( inp_interceptfld_, desc, 0 );
    fillInp( inp_gradientfld_, desc, 1 );
    return true;
}

bool uiAVOPolarAttrib::getOutput( Desc& desc )
{
    fillOutput( desc, outputfld_->getIntValue() );
    return true;
}


