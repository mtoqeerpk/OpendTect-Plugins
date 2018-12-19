/*Copyright (C) 2018 Wayne Mogg All rights reserved.

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
 ________________________________________________________________________

-*/

#include "avopolarattrib.h"
#include "odplugin.h"

mDefODPluginEarlyLoad(AVOPolarAttrib)

mDefODPluginInfo(AVOPolarAttrib)
{
	mDefineStaticLocalObject( PluginInfo, retpi,(
		"AVO Polarization Attribute (base)",
		"AVO Polarization Attribute (base)",
		"Wayne Mogg",
		"6.0",
    	"AVO Polarization, Strength and related attributes for OpendTect v6+",
		PluginInfo::GPL ) );
    return &retpi;
}


mDefODInitPlugin(AVOPolarAttrib)
{
    Attrib::AVOPolarAttrib::initClass();
    return 0;
}

