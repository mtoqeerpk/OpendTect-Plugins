/*Copyright (C) 2018 Wayne Mogg All rights reserved.

This file may be used either under the terms of:

1. The GNU General Public License version 3 or higher, as published by
the Free Software Foundation, or

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef avopolarattrib_h
#define avopolarattrib_h

/*+
________________________________________________________________________

 Author:        Wayne Mogg
 Date:          November 2018
 ________________________________________________________________________

-*/

#include "avopolarattribmod.h"
#include "attribprovider.h"


/*!\brief AVO Polarization Attribute

Calculate AVO polarization, strength and related attributes 

*/


namespace Attrib
{

mClass(AVOPolarAttrib) AVOPolarAttrib : public Provider
{
public:
    static void         initClass();
                        AVOPolarAttrib(Desc&);

    static const char*  attribName()    { return "AVOPolarAttrib"; }
    static const char*	outputStr()     { return "output"; }
    static const char*  gateBGStr()     { return "gateBG"; }
    static const char*  soBGStr()       { return "soBG"; }
    static const char*  gateStr()       { return "gate"; }

    enum OutputType { BackgroundAngle, LocalAngle, AngleDifference, Strength, PolarizationProduct, Quality, BackgroundQuality };

    
protected:
                        ~AVOPolarAttrib();
    static Provider*    createInstance(Desc&);

    bool                allowParallelComputation() const { return false; }

    bool                getInputData(const BinID&,int zintv);
    bool                computeData(const DataHolder&, const BinID& relpos, int z0, int nrsamples, int threadid) const;

    const BinID*        desStepout(int input,int output) const;
    const Interval<int>*    desZSampMargin(int input,int output) const
                            { return &zmargin_; }

    bool                getTrcPos();

    BinID               stepoutBG_;
    Interval<float>     gateBG_;
    Interval<float>     gate_;
    TypeSet<BinID>      trcpos_;
    int                 centertrcidx_;
    int                 outtype_;

    int                 dataidx_;

    ObjectSet<const DataHolder>   interceptdata_;
    ObjectSet<const DataHolder>   gradientdata_;
    int                 interceptdataidx_;
    int                 gradientdataidx_;
};

}; // namespace Attrib


#endif


