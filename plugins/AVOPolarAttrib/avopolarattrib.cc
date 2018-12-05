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
            setOutputValue(output,OutputType::BackgroundAngle, idx, z0, bgAngle[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::LocalAngle)) {
        locAngle.setSize(nrsamples);
        computeLocalAngle(A, B, locAngle);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::LocalAngle, idx, z0, locAngle[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::AngleDifference)) {
        computeAngleDifference();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::AngleDifference, idx, z0, RESULT[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::Strength)) {
        computeStrength();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::Strength, idx, z0, RESULT[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::PolarizationProduct)) {
        computePolarizationProduct();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::PolarizationProduct, idx, z0, RESULT[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::Quality)) {
        computeQuality();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::Quality, idx, z0, RESULT[idx-samplegateBG_.start]);
    }

    if (isOutputEnabled(OutputType::BackgroundQuality)) {
        computeBackgroundQuality();
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::BackgroundQuality, idx, z0, RESULT[idx-samplegateBG_.start]);
    }
    
    return true;
}

void AVOPolarAttrib::computeBackgroundAngle( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& bgAngle )
{
    int ntraces = A.info.getSize(0);
    int sz = A.info.getSize(1);
    Array1DImpl<double> A2(sz);
    A2.setAll(0.0);
    Array1DImpl<double> B2(sz);
    B2.setAll(0.0);
    Array1DImpl<double> AB(sz);
    AB.setAll(0.0);
    
    for (int idx=0; idx<sz; idx++) {
        double A2v = 0.0;
        double ABv = 0.0;
        double B2v = 0.0;
        for (int trcidx=0; trcidx<ntraces; trcidx++) {
            A2v += A.get(trcidx, idx) * A.get(trcidx, idx);
            ABv += A.get(trcidx, idx) * B.get(trcidx, idx);
            B2v += B.get(trcidx, idx) * B.get(trcidx, idx);
        }
        A2[idx] = A2v;
        B2[idx] = B2v;
        AB[idx] = ABv;
    }
    
    int nrsamples = sz - samplegateBG_.width();
    bgAngle.setAll(0.0);
    double A2v = 0.0;
    double ABv = 0.0;
    double B2v = 0.0;
    for (int idx=0; idx<nrsamples; idx++) {
        if (idx==0) {
            for (int i=0; i<samplegateBG_.width(); i++) {
                A2v += A2[i];
                ABv += AB[i];
                B2v += B2[i];
            }
        } else {
            int isub = idx-1;
            int iadd = idx + samplegateBG_.width(); 
            A2v += A2[iadd] - A2[isub];
            B2v += B2[iadd] - B2[isub];
            ABv += AB[iadd] - AB[isub];
        }
        
        double A2mB2 = A2v - B2v;
        double d = sqrt(4.0*ABv*ABv + A2mB2*A2mB2);
        bgAngle[idx] = atan2(2.0*ABv, A2mB2+d);
    }
}

void AVOPolarAttrib::computeLocalAngle( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& locAngle )
{
    int sz = A.info.getSize(1);
    Array1DImpl<double> A2(sz);
    A2.setAll(0.0);
    Array1DImpl<double> B2(sz);
    B2.setAll(0.0);
    Array1DImpl<double> AB(sz);
    AB.setAll(0.0);
    for (int idx=0; idx<sz; idx++) {
        A2[idx] = A.get(centertrcidx_, idx) * A.get(centertrcidx_, idx);
        B2[idx] = B.get(centertrcidx_, idx) * B.get(centertrcidx_, idx);
        AB[idx] = A.get(centertrcidx_, idx) * B.get(centertrcidx_, idx);
    }
    
    int nrsamples = sz - samplegateBG_.width();
    locAngle.setAll(0.0);
    double A2v = 0.0;
    double B2v = 0.0;
    double ABv = 0.0;
    for (int idx=-samplegateBG_.start; idx<nrsamples-samplegateBG_.start; idx++) {
        if (idx == -samplegateBG_.start) {
            for (int i=samplegate_.start+idx; i<samplegate_.stop+idx; i++) {
                A2v += A2[i];
                B2v += A2[i];
                ABv += A2[i];
            }
        } else {
            int isub = idx + samplegate_.start - 1;
            int iadd - idx + samplegate_.stop;
            A2v += A2[iadd] - A2[isub];
            B2v += B2[iadd] - B2[isub];
            ABv += AB[iadd] - AB[isub];
        }
        
        double A2mB2 = A2v - B2v;
        double d = sqrt(4.0*ABv*ABv + A2mB2*A2mB2);
        locAngle[idx] = atan2(2.0*ABv, A2mB2+d);
    }
}

void AVOPolarAttrib::computeStrength( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& strength )
{
    
    int minIdx = 0;
    int maxIdx = 0;
    float minVal = A.get(centertrcidx_, 0);
    float maxVal = min;
    for (int idx=0; idx<sz; idx++) {
        float Aval = A.get(centertrcidx_, idx);
        minIdx = (Aval<minVal)? idx : minIdx;
        maxIdx = (Aval>maxVal)? idx : maxIdx;
        minVal = (Aval<minVal)? Aval : minVal;
        maxVal = (Aval>maxVal)? Aval : maxVal;
    }

    int sz = A.info.getSize(1);
    int nrsamples = sz - samplegateBG_.width();
    strength.setAll(0.0);
   for (int idx=0; idx<nrsamples; idx++) {
        if (idx==0) {
            for (int i=0; i<samplegate_.width(); i++) {
                float Aval = A.get(centertrcidx_, i);
                minIdx = (Aval<minVal)? i : minIdx;
                maxIdx = (Aval>maxVal)? i : maxIdx;
                minVal = (Aval<minVal)? Aval : minVal;
                maxVal = (Aval>maxVal)? Aval : maxVal;
            }
        } else {
            if (minIdx == idx-1) {
                for (int i=0; i<samplegate_.width(); i++) {
                    float Aval = A.get(centertrcidx_, i+idx);
                    minIdx = (Aval<minVal)? i+idx : minIdx;
                    minVal = (Aval<minVal)? Aval : minVal;
                }
            }
            if (maxIdx == idx-1) {
                for (int i=0; i<samplegate_.width(); i++) {
                    float Aval = A.get(centertrcidx_, i+idx);
                    maxIdx = (Aval>maxVal)? i+idx : maxIdx;
                    maxVal = (Aval>maxVal)? Aval : maxVal;
                }
            }
            float Aval = A.get(centertrcidx_, idx+samplegate_.width());
            if (Aval<minVal) {
                minIdx = idx+samplegate_.width();
                minVal = Aval;
            }
            if (Aval>maxVal) {
                maxIdx = idx+samplegate_.width();
                maxVal = Aval;
            }
        }
        float Amin = A.get(minIdx);
        float Bmin = B.get(minIdx);
        float Amax = A.get(maxIdx);
        float Bmax = B.get(maxIdx); 
        strength.set( idx, sqrt(Amin*Amin+Bmin+Bmin)+sqrt(Amax*Amax+Bmax*Bmax));
        }
}

// namespace Attrib

