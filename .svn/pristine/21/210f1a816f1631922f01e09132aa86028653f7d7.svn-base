/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CRYPTCELLCYCLEMODEL
#define CRYPTCELLCYCLEMODEL

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A crypt cell cycle model
 */
class CryptCellCycleModel : public AbstractSimpleGenerationBasedCellCycleModel
{
    friend class TestCryptCellCycleModels;

private:

    /*
     * 1 - Pedegrie Wnt
     * 2 - Spatial Wnt
     * 3 - Spatial Wnt at birth
     * 4 - Mutant
     */
    unsigned mCellProliferationModel;


    /* Aiming for balance of Healthy and Mutant cells if DOUBLE_UNSET then not using this feature */
    double mMutantHealthyRatio;

    /* Helper variable to store the healthy and mutant CCM if using mutant proportion */
    unsigned mHealthyCellProliferationModel;

    unsigned mMutantCellProliferationModel;


    /*
     * If cell is under mQuiescentVolumeFraction then it is held in G1 Phase
     */
    bool mIsContactInhibitionCellCycleDuration;

    bool mIsWntDependentCellCycleDuration;

    /**
     * Non-dimensionalized Wnt threshold, above which cells progress through the cell cycle.
     */
    double mWntThreshold;

    double mMutantWntThreshold;

    /**
     * Store the level of Wnt at the birth of the cell
     */
    double mInitialWntLevel;

    /**
    * The fraction of the cells' equilibrium volume in G1 phase below which these cells are quiescent.
    */
    double mQuiescentVolumeFraction;

    /**
    * The fraction of the mutant cells' equilibrium volume in G1 phase below which these cells are quiescent.
    */
    double mMutantQuiescentVolumeFraction;

    /**
    * The cell equilibrium volume while in G1 phase.
    */
    double mEquilibriumVolume;

    /**
    * The time when the current period of quiescence began.
    */
    double mCurrentQuiescentOnsetTime;

    /**
    * How long the current period of quiescence has lasted.
    * Has units of hours.
    */
    double mCurrentQuiescentDuration;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationBasedCellCycleModel>(*this);
        archive & mCellProliferationModel;
        archive & mMutantHealthyRatio;
		archive & mHealthyCellProliferationModel;
		archive & mMutantCellProliferationModel;
        archive & mIsContactInhibitionCellCycleDuration;
        archive & mIsWntDependentCellCycleDuration;

        archive & mWntThreshold;
        archive & mMutantWntThreshold;
        archive & mInitialWntLevel;

        archive & mQuiescentVolumeFraction;
        archive & mMutantQuiescentVolumeFraction;
        archive & mEquilibriumVolume;
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescentOnsetTime;

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;

    }

protected:

    /**
     * Set the G1 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     *
     * This is overridden here to implement the different CCDs
     */
    void SetG1Duration();

    /**
     * Get the Wnt level experienced by the cell.
     */
    double GetWntLevel();

public:

    /**
     * Constructor. Note, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    CryptCellCycleModel();

    /** Overridden ResetForDivision() method. So doesn't call generation stuff for spatial models*/
    void ResetForDivision();

    /**
     * Overridden UpdateCellCyclePhase() method.
     *
     * This is where we implement different models for cell differentiation.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overriddent InitialiseDaughterCell method to set the CCM based on
     * the target Mutant healthy ratio if being used.
     *
     * Calls parent method to Set the new cell's G1 duration once it has been created after division.
     * The duration will be based on cell type.
     */
    void InitialiseDaughterCell();

    /**
     * @param cellProliferationModel
     */
    void SetCellProliferationModel(unsigned cellProliferationModel);

    /**
     * @return mCellProliferationModel
     */
    unsigned GetCellProliferationModel();

    /**
	 * @param mMutantHealthyRatio
	 */
	void SetMutantHealthyRatio(double mutantHealthyRatio);

    /**
     * @param healthyCellProliferationModel
     */
    void SetHealthyCellProliferationModel(unsigned healthyCellProliferationModel);

    /**
     * @param mutantCellProliferationModel
     */
    void SetMutantCellProliferationModel(unsigned mutantCellProliferationModel);

    /**
     * @param isContactInhibitionCellCycleDuration
     */
    void SetIsContactInhibitionCellCycleDuration(bool isContactInhibitionCellCycleDuration);

    /**
     * @return mIsContactInhibitionCellCycleDuration
     */
    bool GetIsContactInhibitionCellCycleDuration();

    /**
     * @param isWntDependentCellCycleDuration
     */
    void SetIsWntDependentCellCycleDuration(bool isWntDependentCellCycleDuration);

    /**
     * @return mIsWntDependentCellCycleDuration
     */
    bool GetIsWntDependentCellCycleDuration();

    /**
     * @return mWntThreshold if cell is healthy and mMutantWntThreshold if the cells is mutant and were using SetMutantHealthyRatio
     */
    double GetWntThreshold();

    /**
     * Set mWntThreshold.
     *
     * @param wntTransitThreshold the value of mWntThreshold
     */
    void SetWntThreshold(double wntTransitThreshold);

    /**
     * @return mMutantWntThreshold
     */
    double GetMutantWntThreshold();

    /**
     * Set mWntThreshold.
     *
     * @param wntTransitThreshold the value of mWntThreshold
     */
    void SetMutantWntThreshold(double wntMutantTransitThreshold);

    /**
     * Set The initial level of Wnt the cell experiences
     * @param initialWntLevel -the level of wnt
     */
    void SetInitialWntLevel(double initialWntLevel);

    /**
     * @param quiescentVolumeFraction
     */
    void SetQuiescentVolumeFraction(double quiescentVolumeFraction);

    /**
     * @return mQuiescentVolumeFraction
     */
    double GetQuiescentVolumeFraction();

    /**
     * @param mutantQuiescentVolumeFraction
     */
    void SetMutantQuiescentVolumeFraction(double mutantQuiescentVolumeFraction);

    /**
     * @return mMutantQuiescentVolumeFraction
     */
    double GetMutantQuiescentVolumeFraction();

    /**
     * @param equilibriumVolume
     */
    void SetEquilibriumVolume(double equilibriumVolume);

    /**
     * @return mEquilibriumVolume
     */
    double GetEquilibriumVolume();

    /**
    * Set method for mCurrentQuiescentDuration.
    *
    * @param currentQuiescentDuration the new value of mCurrentQuiescentDuration
    */
    void SetCurrentQuiescentDuration(double currentQuiescentDuration);

    /**
     * @return mCurrentQuiescentDuration
     */
    double GetCurrentQuiescentDuration();

    /**
    * Set method for mCurrentQuiescentOnsetTime.
    *
    * @param currentQuiescentOnsetTime the new value of mCurrentQuiescentOnsetTime
    */
    void SetCurrentQuiescentOnsetTime(double currentQuiescentOnsetTime);

    /**
     * @return mCurrentQuiescentOnsetTime
     */
    double GetCurrentQuiescentOnsetTime();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CryptCellCycleModel)

#endif /*CRYPTCELLCYCLEMODEL*/
