/*Copyright (C) 2014 Wayne Mogg. All rights reserved.
 * 
 T hi*s file may be used either under the terms of:
 
 1. The GNU General Public License version 3 or higher, as published by
 the Free Software Foundation, or
 
 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*+
 _ __*_____________________________________________________________________
 
 Author:        Wayne Mogg
 Date:          January 2014
 ________________________________________________________________________
 
 -*/

#include "uimenu.h"
#include "uiodmain.h"
#include "odplugin.h"

#include "uiavoattrib.h"


mDefODPluginInfo(uiAVOAttrib)
{
	mDefineStaticLocalObject( PluginInfo, retpi,(
		"AVO Attribute (UI)",
		"AVO Attribute (UI)",
		"Wayne Mogg",
		"6.0",
		"",
		PluginInfo::GPL ) );	
    return &retpi;
}


mDefODInitPlugin(uiAVOAttrib)
{
	uiAVOAttrib::initClass();
	return 0;
}
