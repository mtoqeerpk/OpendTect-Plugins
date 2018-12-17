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
#include "attribdescset.h"
#include "attribfactory.h"
#include "attribparam.h"
#include "arrayndimpl.h"

#include "windowedOps.h"


namespace Attrib
{

mAttrDefCreateInstance(AVOPolarAttrib)    
    
void AVOPolarAttrib::initClass()
{
	mAttrStartInitClass

//    EnumParam* output = new EnumParam( outputStr() );
//    output->addEnum( "Background Angle" );
//    output->addEnum( "Local Angle" );
//    output->addEnum( "Angle Difference" );
//    output->addEnum( "Strength" );
//    output->addEnum( "Polarization Product" );
//    output->addEnum( "Quality" );
//    desc->addParam( output );
    
    ZGateParam* gateBG = new ZGateParam( gateBGStr() );
    gateBG->setLimits( Interval<float>(-1000,1000) );
    gateBG->setDefaultValue( Interval<float>(-100, 100) );
    desc->addParam( gateBG );

    BinIDParam* stepout = new BinIDParam( soBGStr() );
    stepout->setDefaultValue( BinID(1,1) );
    desc->addParam( stepout );
    
    desc->addInput( InputSpec("Intercept cube",true) );
    desc->addInput( InputSpec("Gradient cube",true) );

    desc->setNrOutputs( Seis::UnknowData, 6 );
    mAttrEndInitClass
}

AVOPolarAttrib::AVOPolarAttrib( Desc& desc )
    : Provider( desc )
{
    if ( !isOK() ) return;

//    mGetEnum( output_, outputStr() );
    mGetFloatInterval( gateBG_, gateBGStr() );
    gateBG_.scale( 1.f/zFactor() );
    mGetFloatInterval( gate_, gateStr() );
    gate_.scale( 1.f/zFactor() );
    mGetBinID( stepoutBG_, soBGStr() );

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
    for ( bid.inl()=-stepoutBG_.inl(); bid.inl()<=stepoutBG_.inl(); bid.inl()++ )
    {
        for ( bid.crl()=-stepoutBG_.crl(); bid.crl()<=stepoutBG_.crl(); bid.crl()++ )
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

    for ( int idx=0; idx<trcpos_.size(); idx++ )
    {
        const DataHolder* data = inputs_[1]->getData( relpos+trcpos_[idx]*bidstep, zintv );
        if ( !data ) {
            const BinID pos = relpos + trcpos_[centertrcidx_]*bidstep;
            data = inputs_[1]->getData( pos, zintv );
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
    
    for (int trcidx=0; trcidx<ntraces; trcidx++) {
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
    Array1DImpl<float> correlationCoef(0);
    Array1DImpl<float> angleDiff(0);
    Array1DImpl<float> strength(0);
    
    if (isOutputEnabled(OutputType::BackgroundAngle)) {
        computeBackgroundAngle(A, B, bgAngle);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::BackgroundAngle, idx, z0, bgAngle[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::LocalAngle)) {
        computeLocalAngle(A, B, locAngle, correlationCoef);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::LocalAngle, idx, z0, locAngle[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::AngleDifference)) {
        if (bgAngle.info().getSize(0) == 0)
            computeBackgroundAngle(A, B, bgAngle);
        if (locAngle.info().getSize(0) == 0)
            computeLocalAngle(A, B, locAngle, correlationCoef);
        computeAngleDifference(bgAngle, locAngle, angleDiff);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::AngleDifference, idx, z0, angleDiff[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::Strength)) {
        computeStrength( A, B, strength);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::Strength, idx, z0, strength[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::PolarizationProduct)) {
        if (strength.info().getSize(0) == 0)
            computeStrength( A, B, strength );
        if (angleDiff.info().getSize(0) == 0) {
            if (bgAngle.info().getSize(0) == 0)
                computeBackgroundAngle(A, B, bgAngle);
            if (locAngle.info().getSize(0) == 0)
                computeLocalAngle(A, B, locAngle, correlationCoef);
            computeAngleDifference(bgAngle, locAngle, angleDiff);
        }
        Array1DImpl<float> polProd(0);
        computePolarizationProduct(angleDiff, strength, polProd);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::PolarizationProduct, idx, z0, polProd[idx-samplegateBG_.start]);
    }
    
    if (isOutputEnabled(OutputType::Quality)) {
        if (correlationCoef.info().getSize(0) == 0)
            computeLocalAngle(A, B, locAngle, correlationCoef);
        for (int idx=0; idx<nrsamples; idx++)
            setOutputValue(output,OutputType::Quality, idx, z0, correlationCoef[idx-samplegateBG_.start]);
    }

    
    return true;
}

void AVOPolarAttrib::computeBackgroundAngle( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& bgAngle ) const
{
    int ntraces = A.info().getSize(0);
    int sz = A.info().getSize(1);
    bgAngle.setSize(sz);
    
    Array1DImpl<float> A2(sz);
    A2.setAll(0.0);
    Array1DImpl<float> B2(sz);
    B2.setAll(0.0);
    Array1DImpl<float> AB(sz);
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
        A2.set(idx, A2v);
        B2.set(idx, B2v);
        AB.set(idx, ABv);
    }
    
    Array1DImpl<float> A2win(sz);
    windowedOps::sum( A2, samplegateBG_.width(), A2win );
    Array1DImpl<float> B2win(sz);
    windowedOps::sum( B2, samplegateBG_.width(), B2win );
    Array1DImpl<float> ABwin(sz);
    windowedOps::sum( AB, samplegateBG_.width(), ABwin );
    
    for (int idx=0; idx<sz; idx++) {
        double ABv = ABwin.get(idx);
        double A2mB2 = A2win.get(idx) - B2win.get(idx);
        double d = sqrt(4.0*ABv*ABv + A2mB2*A2mB2);
        bgAngle.set(idx, atan2(2.0*ABv, A2mB2+d));
    }
}
        
void AVOPolarAttrib::computeLocalAngle( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& locAngle, Array1DImpl<float>& coeff ) const
{
    int sz = A.info().getSize(1);
    locAngle.setSize(sz);
    coeff.setSize(sz);

    Array1DImpl<float> A2(sz);
    A2.setAll(0.0);
    Array1DImpl<float> B2(sz);
    B2.setAll(0.0);
    Array1DImpl<float> AB(sz);
    AB.setAll(0.0);
    for (int idx=0; idx<sz; idx++) {
        A2.set(idx, A.get(centertrcidx_, idx) * A.get(centertrcidx_, idx));
        B2.set(idx, B.get(centertrcidx_, idx) * B.get(centertrcidx_, idx));
        AB.set(idx, A.get(centertrcidx_, idx) * B.get(centertrcidx_, idx));
    }
    
    Array1DImpl<float> A2win(sz);
    windowedOps::sum( A2, samplegate_.width(), A2win );
    Array1DImpl<float> B2win(sz);
    windowedOps::sum( B2, samplegate_.width(), B2win );
    Array1DImpl<float> ABwin(sz);
    windowedOps::sum( AB, samplegate_.width(), ABwin );
    
    for (int idx=0; idx<sz; idx++) {
        double ABv = ABwin.get(idx);
        double A2mB2 = A2win.get(idx) - B2win.get(idx);
        double d = sqrt(4.0*ABv*ABv + A2mB2*A2mB2);
        locAngle.set(idx, atan2(2.0*ABv, A2mB2+d));
        coeff.set(idx, ABv / (sqrt(A2win.get(idx)) * sqrt(B2win.get(idx))));
    }
}

void AVOPolarAttrib::computeAngleDifference( const Array1DImpl<float>& bgAngle, const Array1DImpl<float>& locAngle, Array1DImpl<float>& diff ) const
{
    int sz = bgAngle.info().getSize(0);
    diff.setSize(sz);
    for (int idx=0; idx<sz; idx++)
        diff.set(idx, locAngle.get(idx)-bgAngle.get(idx));
}

void AVOPolarAttrib::computeStrength( const Array2DImpl<float>& A, const Array2DImpl<float>& B, Array1DImpl<float>& strength ) const
{
    int sz = A.info().getSize(1);
    strength.setSize(sz);

    Array1DImpl<float> Ac(sz);
    Array1DImpl<float> Bc(sz);
    for (int idx=0; idx<sz; idx++) {
        Ac.set(idx, A.get(centertrcidx_, idx));
        Bc.set(idx, B.get(centertrcidx_, idx));
    }
    
    Array1DImpl<int> Amin(sz);
    windowedOps::minIdx( Ac, samplegate_.width(), Amin );
    Array1DImpl<int> Amax(sz);
    windowedOps::maxIdx( Ac, samplegate_.width(), Amax );
    
    for (int idx=0; idx<sz; idx++) {
        int imin = Amin.get(idx);
        int imax = Amin.get(idx);
        double Amn = Ac.get(imin);
        double Bmn = Bc.get(imin);
        double Amx = Ac.get(imax);
        double Bmx = Bc.get(imax);
        strength.set(idx, sqrt(Amn*Amn+Bmn*Bmn) + sqrt(Amx*Amx+Bmx*Bmx));
    }    
}

void AVOPolarAttrib::computePolarizationProduct( const Array1DImpl<float>& angleDiff, const Array1DImpl<float>& strength, Array1DImpl<float>& polProd ) const
{
    int sz = angleDiff.info().getSize(0);
    polProd.setSize(sz);

    for (int idx=0; idx<sz; idx++) {
        polProd.set(idx, angleDiff.get(idx)*strength.get(idx));
    }    
}
}

// namespace Attrib

