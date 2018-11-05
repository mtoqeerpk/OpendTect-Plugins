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
 ________________________________________________________________________

-*/

#include "uimenu.h"
#include "uiodmain.h"
#include "odplugin.h"

#include "uiavopolarattrib.h"


mDefODPluginInfo(uiAVOPolarAttrib)
{
	mDefineStaticLocalObject( PluginInfo, retpi,(
		"AVO Polarization, Strength and related attributes (UI)",
		"AVO Polarization, Strength and related attributes (UI)",
		"Wayne Mogg",
		"6.0",
		"",
		PluginInfo::GPL ) );
    return &retpi;
}


mDefODInitPlugin(uiAVOPolarAttrib)
{
    uiAVOPolarAttrib::initClass();
   return 0;
}

