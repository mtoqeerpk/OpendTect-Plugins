/*Copyright (C) 2018 Wayne Mogg All rights reserved.
 * 
 * This file may be used either under the terms of:
 * 
 * 1. The GNU General Public License version 3 or higher, as published by
 * the Free Software Foundation, or
 * 
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */
/*+
 * ________________________________________________________________________
 * 
 * Author:        Wayne Mogg
 * Date:          November 2018
 * ________________________________________________________________________
 * 
 * -*/


#include "avopolarattrib.h"
#include "attribdataholder.h"
#include "attribdesc.h"
#include "attribfactory.h"
#include "attribparam.h"
#include "arrayndimpl.h"

#include "math2.h"


namespace Attrib
{

mAttrDefCreateInstance(AVOPolarAttrib)    
    
void AVOPolarAttrib::initClass()
{
	mAttrStartInitClassWithUpdate

    EnumParam* output = new EnumParam( outputStr() );
    output->addEnum( "Background Angle" );
    output->addEnum( "Local Angle" );
    output->addEnum( "Angle Difference" );
    output->addEnum( "Strength" );
    output->addEnum( "Polarization Product" );
    output->addEnum( "Quality" );
    output->addEnum( "Background Quality" );
    desc->addParam( output );
    
    ZGateParam* gateBG = new ZGateParam( gateBGStr() );
    gateBG->setLimits( Interval<float>(-1000,1000) );
    gateBG->setDefaultValue( Interval<float>(-100, 100) );
    desc->addParam( gateBG );

    BinIDParam* stepout = new BinIDParam( soBGStr() );
    stepout->setDefaultValue( BinID(1,1) );
    desc->addParam( stepout );
    
    desc->addOutputDataType( Seis::UnknowData );
    desc->addInput( InputSpec("Intercept cube",true) );
    desc->addInput( InputSpec("Gradient cube",true) );

    mAttrEndInitClass
}

void AVOPolarAttrib::updateDesc( Desc& desc )
{
    BufferString output = desc.getValParam( outputStr() )->getStringValue();
}

AVOPolarAttrib::AVOPolarAttrib( Desc& desc )
    : Provider( desc )
{
    if ( !isOK() ) return;

    mGetEnum( output_, outputStr() );
    mGetFloatInterval( gateBG_, gateBGStr() );
    gateBG_.scale( 1.f/zFactor() );
    mGetFloatInterval( gate_, gateStr() );
    gate_.scale( 1.f/zFactor() );
    mGetBinID( stepoutBG_, soBGStr() )

    samplegateBG_ = Interval<int>( mNINT32(gateBG_.start/refstep_), mNINT32(gateBG_.stop/refstep_) );
    samplegate_ = Interval<int>( mNINT32(gate_.start/refstep_), mNINT32(gate_.stop/refstep_) );

    getTrcPos();
}

bool AVOPolarAttrib::getTrcPos()
{
    trcpos_.erase();
    BinID bid;
    int trcidx = 0;
    centertrcidx_ = 0;
    for ( bid.inl()=-stepout_.inl(); bid.inl()<=stepout_.inl(); bid.inl()++ )
    {
        for ( bid.crl()=-stepout_.crl(); bid.crl()<=stepout_.crl(); bid.crl()++ )
        {
            if ( !bid.inl() && !bid.crl() )
                centertrcidx_ = trcidx;
            trcpos_ += bid;
            trcidx++;
        }
    }
    return true;
}

bool AVOPolarAttrib::getInputData( const BinID& relpos, int zintv )
{
    while ( interceptdata_.size() < trcpos_.size() )
        interceptdata_ += 0;

    while ( gradientdata_.size() < trcpos_.size() )
        gradientdata_ += 0;

    const BinID bidstep = inputs_[0]->getStepoutStep();
    for ( int idx=0; idx<trcpos_.size(); idx++ )
    {
        const DataHolder* data = inputs_[0]->getData( relpos+trcpos_[idx]*bidstep, zintv );
        if ( !data ) {
            const BinID pos = relpos + trcpos_[centertrcidx_]*bidstep;
            data = inputs_[0]->getData( pos, zintv );
            if ( !data ) return false;
        }
        interceptdata_.replace( idx, data );
    }

    bidstep = inputs_[1]->getStepoutStep();
    for ( int idx=0; idx<trcpos_.size(); idx++ )
    {
        const DataHolder* data = inputs_[1]->getData( relpos+trcpos_[idx]*bidstep, zintv );
        if ( !data ) {
            const BinID pos = relpos + trcpos_[centertrcidx_]*bidstep;
            data = inputs_[0]->getData( pos, zintv );
            if ( !data ) return false;
        }
        gradientdata_.replace( idx, data );
    }

    interceptdataidx_ = getDataIndex( 0 );
    gradientdataidx_ = getDataIndex( 1 );

    return true;
}


bool AVOPolarAttrib::computeData( const DataHolder& output, const BinID& relpos, int z0, int nrsamples, int threadid ) const
{
    if ( interceptdata_.isEmpty() || gradientdata_.isEmpty() || output.isEmpty() )
        return false;

    const int ntraces = trcpos_.size();
    const int sz = samplegateBG_.width() + nrsamples;
    
    Array2DImpl<float> A(ntraces, sz), B(ntraces, sz);
    
    for (int trcidx=0; trcidx<ntraces; trc++) {
        const DataHolder* dataA = interceptdata_[trcidx];
        const DataHolder* dataB = gradientdata_[trcidx];
        for (int idx=0; idx<sz; idx++) {
            float val = getInputValue(*dataA, interceptdataidx_, samplegateBG_.start+idx, z0);
            A.set(trcidx, idx, mIsUdf(val)?0.0f:val);
            val = getInputValue(*dataB, gradientdataidx_, samplegateBG_.start+idx, z0);
            B.set(trcidx, idx, mIsUdf(val)?0.0f:val);
        }
    }

    Array1DImpl<float> bgAngle(0);
    Array1DImpl<float> locAngle(0);
    
    if (isOutputEnabled(OutputType::BackgroundAngle)) {
        bgAngle.setSize(nrsamples);
        computeBackgroundAngle(A, B, bgAngle);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::BackgroundAngle, idx, z0, bgAngle[idx]);
    }
    
    if (isOutputEnabled(OutputType::LocalAngle)) {
        locAngle.setSize(nrsamples);
        computeLocalAngle(A, B, locAngle);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::LocalAngle, idx, z0, locAngle[idx]);
    }
    
    if (isOutputEnabled(OutputType::AngleDifference)) {
        computeAngleDifference();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::AngleDifference, idx, z0, RESULT[idx]);
    }
    
    if (isOutputEnabled(OutputType::Strength)) {
        computeStrength();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::Strength, idx, z0, RESULT[idx]);
    }
    
    if (isOutputEnabled(OutputType::PolarizationProduct)) {
        computePolarizationProduct();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::PolarizationProduct, idx, z0, RESULT[idx]);
    }
    
    if (isOutputEnabled(OutputType::Quality)) {
        computeQuality();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::Quality, idx, z0, RESULT[idx]);
    }

    if (isOutputEnabled(OutputType::BackgroundQuality)) {
        computeBackgroundQuality();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::BackgroundQuality, idx, z0, RESULT[idx]);
    }
    
    return true;
}

void AVOPolarAttrib::computeBackgroundAngle( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& bgAngle )
{
    int ntraces = A.info.getSize(0);
    int sz = A.info.getSize(1);
    Array1DImpl<double> A2(sz);
    A2.setAll(0.0);
    
    for (int idx=0; idx<sz; idx++) {
        double A2v = 0.0
        for (int trcidx=0; trcidx<ntraces; trcidx++) {
            A2v += A.get(trcidx, idx);
        }
        A2[idx] = A2v;
    }
}
} // namespace Attrib

