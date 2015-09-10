#include <iostream>
#include "chromoconditions.h"

namespace BioLCCC
{
ChromoConditionsException::ChromoConditionsException(std::string message):
        BioLCCCException(message) {};

ChromoConditions::ChromoConditions(double iColumnLength,
                                   double iColumnDiameter,
                                   double iColumnPoreSize,
                                   Gradient iGradient,
                                   double iSecondSolventConcentrationA,
                                   double iSecondSolventConcentrationB,
                                   double iDelayTime,
                                   double iFlowRate,
                                   double iDV,
                                   double iColumnRelativeStrength,
                                   double iColumnVpToVtot,
                                   double iColumnPorosity,
                                   double iTemperature)
                                   throw(ChromoConditionsException)
{
    // Set an empty gradient to prevent recalculation of SSConcentrations.
    mGradient = Gradient();
    setMixingCorrection(false);
    setColumnLength(iColumnLength);
    setColumnDiameter(iColumnDiameter);
    setColumnPoreSize(iColumnPoreSize);
    setColumnVpToVtot(iColumnVpToVtot);
    setColumnPorosity(iColumnPorosity);
    setTemperature(iTemperature);
    setColumnRelativeStrength(iColumnRelativeStrength);
    setFlowRate(iFlowRate);
    setDV(iDV);
    setDelayTime(iDelayTime);
    setSecondSolventConcentrationA(iSecondSolventConcentrationA);
    setSecondSolventConcentrationB(iSecondSolventConcentrationB);
    setGradient(iGradient);
}

double ChromoConditions::columnLength() const
{
    return mColumnLength;
}

void ChromoConditions::setColumnLength(double newColumnLength)
    throw(ChromoConditionsException)
{
    if (newColumnLength < 0.0)
    {
        throw(ChromoConditionsException("The new column length is negative."));
    }
    mColumnLength = newColumnLength;
    recalculateVolumes();
    recalculateSSConcentrations();
}

double ChromoConditions::columnDiameter() const
{
    return mColumnDiameter;
}

void ChromoConditions::setColumnDiameter(double newColumnDiameter)
    throw(ChromoConditionsException)
{
    if (newColumnDiameter < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column diameter is negative."));
    }
    mColumnDiameter = newColumnDiameter;
    recalculateVolumes();
    recalculateSSConcentrations();
}

double ChromoConditions::columnPoreSize() const
{
    return mColumnPoreSize;
}

void ChromoConditions::setColumnPoreSize(double newColumnPoreSize)
    throw(ChromoConditionsException)
{
    if (newColumnPoreSize < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column pore size is negative."));
    }
    mColumnPoreSize = newColumnPoreSize;
}

double ChromoConditions::columnVpToVtot() const
{
    return mColumnVpToVtot;
}

void ChromoConditions::setColumnVpToVtot(double newColumnVpToVtot)
    throw(ChromoConditionsException)
{
    if (newColumnVpToVtot < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column VpToVtot is negative."));
    }
    if (newColumnVpToVtot > 1.0)
    {
        throw(ChromoConditionsException(
            "The new column VpToVtot is greater than 1.0."));
    }
    mColumnVpToVtot = newColumnVpToVtot;
    recalculateVolumes();
    recalculateSSConcentrations();
}

double ChromoConditions::columnPorosity() const
{
    return mColumnPorosity;
}

void ChromoConditions::setColumnPorosity(double newColumnPorosity)
    throw(ChromoConditionsException)
{
    if (newColumnPorosity < 0.0)
    {
        throw(ChromoConditionsException(
            "The new column porosity is negative."));
    }
    if (newColumnPorosity > 1.0)
    {
        throw(ChromoConditionsException(
            "The new column porosity is greater than 1.0"));
    }
    mColumnPorosity = newColumnPorosity;
    recalculateVolumes();
    recalculateSSConcentrations();
}

double ChromoConditions::columnTotalVolume() const
{
    return mColumnTotalVolume;
}

double ChromoConditions::columnInterstitialVolume() const
{
    return mColumnInterstitialVolume;
}

double ChromoConditions::columnPoreVolume() const
{
    return mColumnPoreVolume;
}

double ChromoConditions::temperature() const
{
    return mTemperature;
}

void ChromoConditions::setTemperature(double newTemperature)
    throw(ChromoConditionsException)
{
    if (newTemperature < 0.0)
    {
        throw(ChromoConditionsException("The new temperature is negative."));
    }
    mTemperature = newTemperature;
}

double ChromoConditions::columnRelativeStrength() const
{
    return mColumnRelativeStrength;
}

void ChromoConditions::setColumnRelativeStrength(
    double newColumnRelativeStrength)
{
    mColumnRelativeStrength = newColumnRelativeStrength;
}

double ChromoConditions::flowRate() const
{
    return mFlowRate;
}

void ChromoConditions::setFlowRate(double newFlowRate)
    throw(ChromoConditionsException)
{
    if (newFlowRate < 0.0)
    {
        throw(ChromoConditionsException("The new flow rate is negative."));
    }
    mFlowRate = newFlowRate;
    recalculateSSConcentrations();
}

double ChromoConditions::dV() const
{
    if (mDV == 0.0)
    {
        return flowRate() / 20.0;
    }
    else
    {
        return mDV;
    }
}

void ChromoConditions::setDV(double newDV)
    throw(ChromoConditionsException)
{
    if (newDV < 0.0)
    {
        throw(ChromoConditionsException("The new dV is negative."));
    }
    mDV = newDV;
    recalculateSSConcentrations();
}

double ChromoConditions::delayTime() const
{
    return mDelayTime;
}

void ChromoConditions::setDelayTime(double newDelayTime)
{
    mDelayTime = newDelayTime;
}

double ChromoConditions::secondSolventConcentrationA() const
{
    return mSecondSolventConcentrationA;
}

void ChromoConditions::setSecondSolventConcentrationA(
    double newSecondSolventConcentrationA)
    throw(ChromoConditionsException)
{
    if (newSecondSolventConcentrationA < 0.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is negative."));
    }
    if (newSecondSolventConcentrationA > 100.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is greater than 100%."));
    }
    mSecondSolventConcentrationA = newSecondSolventConcentrationA;
}

double ChromoConditions::secondSolventConcentrationB() const
{
    return mSecondSolventConcentrationB;
}

void ChromoConditions::setSecondSolventConcentrationB(
    double newSecondSolventConcentrationB)
    throw(ChromoConditionsException)
{
    if (newSecondSolventConcentrationB < 0.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is negative."));
    }
    if (newSecondSolventConcentrationB > 100.0)
    {
        throw(ChromoConditionsException(
            "The new concentration of second solvent in the component A "
            "is greater than 100%."));
    }
    mSecondSolventConcentrationB = newSecondSolventConcentrationB;
}

Gradient ChromoConditions::gradient() const
{
    return mGradient;
}

void ChromoConditions::setGradient(Gradient newGradient)
    throw(ChromoConditionsException)
{
    if (newGradient.size() < 2)
    {
        throw ChromoConditionsException(
            "The gradient must contain at least two points.");
    }
    mGradient = newGradient;
    recalculateSSConcentrations();
}

const std::vector<double> & ChromoConditions::SSConcentrations() const
{
    return mSSConcentrations;
}

void ChromoConditions::recalculateVolumes()
{
    mColumnTotalVolume = columnDiameter() * columnDiameter() * 3.1415 / 4.0 *
                         columnLength() / 1000.0;

    mColumnInterstitialVolume = mColumnTotalVolume * 
                                (columnPorosity() -  columnVpToVtot());

    mColumnPoreVolume = mColumnTotalVolume * columnVpToVtot();
}

bool ChromoConditions::mixingCorrection() const
{
    return mMixingCorrection;
}
    
void ChromoConditions::setMixingCorrection(bool flag)
{
    mMixingCorrection = flag;
    recalculateSSConcentrations();
}

void ChromoConditions::recalculateSSConcentrations()
{
    mSSConcentrations.clear();
    if (gradient().empty())
    {
        return;
    }
    
    // If the gradient is isocratic, add a single point only.
    if ((gradient().size() == 2) 
         && (gradient().front().concentrationB() 
             == gradient().back().concentrationB()))
    {
        mSSConcentrations.push_back(
            (100.0 - gradient().front().concentrationB()) / 100.0 
            * secondSolventConcentrationA() 
            + gradient().front().concentrationB() / 100.0
            * secondSolventConcentrationB());
        return;
    }

    double secondSolventConcentrationPump = 0.0;
    double time = dV() / 2.0 / flowRate();
    int segmentNum = -1;
    double localSlope = 0.0;
    double initialSSConcentration, finalSSConcentration, initialTime,
           finalTime, pumpedConcentration;

    while (true)
    {
        if (time > gradient()[segmentNum+1].time())
        {
            if (segmentNum < int(gradient().size() - 1))
            {
                segmentNum += 1;
                initialSSConcentration =
                    (100.0 - gradient()[segmentNum].concentrationB()) / 100.0 
                    * secondSolventConcentrationA() 
                    + gradient()[segmentNum].concentrationB() / 100.0
                    * secondSolventConcentrationB();
                finalSSConcentration =
                    (100.0 - gradient()[segmentNum+1].concentrationB()) / 100.0 
                    * secondSolventConcentrationA() 
                    + gradient()[segmentNum+1].concentrationB() / 100.0
                    * secondSolventConcentrationB();
                initialTime = gradient()[segmentNum].time();
                finalTime = gradient()[segmentNum+1].time();
                localSlope = (finalSSConcentration - initialSSConcentration) 
                    / (finalTime - initialTime);
            }
            else
            {
                break;
            }
        }

        pumpedConcentration = localSlope * (time - initialTime) 
                              + initialSSConcentration;

        // If mixingCorrection is enabled calculate second solvent concentrations
        // from the equation 
        // d[SS] / dt = flowRate / (V0 + VP) * ([SS]pump - [SS])
        // Otherwise, [SS] == [SS]pump
        if (mixingCorrection())
        {
            // Special case for the first time point.
            if (time == dV() / 2.0 / flowRate())
            {
                mSSConcentrations.push_back(
                    initialSSConcentration + 
                    + dV() / (columnInterstitialVolume() + columnPoreVolume())
                    * (pumpedConcentration - mSSConcentrations.back()) / 2.0);
            }
            else
            {
                mSSConcentrations.push_back(
                    mSSConcentrations.back() 
                    + dV() / (columnInterstitialVolume() + columnPoreVolume())
                    * (pumpedConcentration - mSSConcentrations.back()));
            }
        }
        else
        {
            mSSConcentrations.push_back(pumpedConcentration);
        }
        time += dV() / flowRate();
    }
}
}

