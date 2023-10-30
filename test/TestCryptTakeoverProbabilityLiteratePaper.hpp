#ifndef TESTCRYPTTAKEOVERPROBABILITYLITERATEPAPER_HPP_
#define TESTCRYPTTAKEOVERPROBABILITYLITERATEPAPER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "OffLatticeSimulationWithMonoclonalStoppingEvent.hpp"
#include "VolumeTrackingModifier.hpp"

#include "CryptCellsGenerator.hpp"
#include "CellsGenerator.hpp"

#include "CylindricalHoneycombMeshGenerator.hpp"

#include "CryptCellCycleModel.hpp"

#include "PanethCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellAncestorWriter.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionSpringForce.hpp"
#include "RepulsionForce.hpp"
#include "CellRetainerForce.hpp"

#include "CryptSimulationBoundaryCondition.hpp"
#include "CryptGeometryBoundaryCondition3d.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"


#include "SloughingCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCryptTakeoverProbabilityLiteratePaper : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void Test3DCrypt() throw (Exception)
    {

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-run_index"));
        unsigned start_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-run_index");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_runs"));
        unsigned num_runs = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-num_runs");


    	unsigned num_sweeps= 11;
		double percent_mutant_array[11] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

        double time_to_steady_state = 100.0;
        double time_after_mutations = 10000.0;

    	// Optimal Healthy Model
    	unsigned healthy_cell_proliferation_model = 3u; // Spatial Wnt at birth
    	bool healthy_wnt_dependend_ccd = true;
		double healthy_wnt_thresh = 0.6;
		double healthy_CI = 0.9;

		//Optimal Mutant model, as above with
		unsigned mutant_cell_proliferation_model = 3u; // Spatial Wnt at birth
		bool mutant_wnt_dependend_ccd = true;
	    double mutant_wnt_thresh = 0.5;
		double mutant_CI = 0.6;

        // Crypt Setup
        double cell_radius = 3.5;//3.5;
        double crypt_length = 70; //70
        double crypt_radius = 8.0/M_PI*6.0; // Choosing same dimensions as for halted migration paper
        // For this size domain there are about 75 cells in the bottom hemisphere so this makes about 20% Paneth cells

		unsigned num_paneth_cells = 15;//15;
		unsigned num_stem_cells = 60; // 60;
		unsigned num_cells = num_paneth_cells + num_stem_cells;

		double stem_retainer_force_magnitude = 7.5*10;
		double paneth_retainer_force_magnitude = 7.5*10;

		// Loop over the random seed.
		for(unsigned sim_index=start_index; sim_index < start_index + num_runs; sim_index++)
		{
			std::cout << " Run number " << sim_index << "... " << std::flush;

			// Reseed the random number generator
			RandomNumberGenerator::Instance()->Reseed(100*sim_index);

			for (unsigned param_index = 0; param_index<num_sweeps;   param_index ++)
			{
				double percent_mutant = percent_mutant_array[param_index];

				PRINT_2_VARIABLES(param_index,
									percent_mutant);


				// Create some starter nodes
				std::vector<Node<3>*> nodes;
				for(unsigned node_index= 0;  node_index<num_cells; node_index++)
				{
					double x = crypt_radius/2.0 * sin(node_index*2.0*M_PI/num_cells);
					double y = crypt_radius/2.0 * cos(node_index*2.0*M_PI/num_cells);
					double z = 0.0;
					nodes.push_back(new Node<3>(node_index, false, x, y, z));
				}

				// Convert this to a NodesOnlyMesh
				NodesOnlyMesh<3> mesh;
				mesh.ConstructNodesWithoutMesh(nodes,cell_radius*3.0);


				// Create cells
				std::vector<CellPtr> cells;

				CellsGenerator<CryptCellCycleModel, 3> cells_generator;
				cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

				//Change properties of the ccm
				MAKE_PTR(PanethCellProliferativeType, p_paneth_type);
				boost::shared_ptr<AbstractCellProperty> p_mutant_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
				boost::shared_ptr<AbstractCellProperty> p_paneth_mutant_state(CellPropertyRegistry::Instance()->Get<ApcOneHitCellMutationState>());

				for (unsigned cell_index= 0;  cell_index<cells.size(); cell_index++)
				{
					cells[cell_index]->GetCellData()->SetItem("Radius", cell_radius);

					// Specify CCM
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetCellProliferationModel(healthy_cell_proliferation_model);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetIsContactInhibitionCellCycleDuration(true);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetIsWntDependentCellCycleDuration(healthy_wnt_dependend_ccd);

                    // Set some default CCD parameters So total CCM is U[10,14] and (U[22,26] at base if variable)
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetMDuration(4.0);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetSDuration(4.0);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetG2Duration(2.0);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetTransitCellG1Duration(2.0);  // so total CCM is U[10,14] at threshold
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetStemCellG1Duration(14.0);  // so total CCM is U[10,14] at base

                    // All cells are intially healthy cells
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetWntThreshold(healthy_wnt_thresh);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetMutantWntThreshold(healthy_wnt_thresh);

                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetQuiescentVolumeFraction(healthy_CI);

                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetMaxTransitGenerations(UINT_MAX);
                    dynamic_cast<CryptCellCycleModel*>(cells[cell_index]->GetCellCycleModel())->SetEquilibriumVolume(M_PI*4.0/3.0*cell_radius*cell_radius*cell_radius);
				}

				// Make some cells paneth cells
				for(unsigned cell_index= 0;  cell_index<num_paneth_cells; cell_index++)
				{
						unsigned temp_index = cell_index * num_cells/num_paneth_cells;
						if (temp_index > num_cells)
						{
							temp_index = num_cells;
						}
						cells[temp_index]->SetCellProliferativeType(p_paneth_type);
						// Give them a mutation to make them easier to track (doesn't actually do anything but label the cells)
						cells[temp_index]->SetMutationState(p_paneth_mutant_state);
				}

				// Create cell population
				NodeBasedCellPopulation<3> cell_population(mesh, cells);
				cell_population.SetUseVariableRadii(true);

				// Output data
				cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
				cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
				cell_population.AddCellWriter<CellAgesWriter>();
				cell_population.AddCellWriter<CellVolumesWriter>();
				cell_population.AddCellWriter<CellProliferativeTypesWriter>();
				cell_population.AddCellWriter<CellMutationStatesWriter>();
				cell_population.AddPopulationWriter<NodeVelocityWriter>();
				cell_population.AddCellWriter<CellAncestorWriter>();

				cell_population.SetAbsoluteMovementThreshold(50.0);


				// Create an instance of a Wnt concentration NOTE DO THIS BEFORE THE SIMULATION OTHERWISE CELLS CANT INITIALISE
				WntConcentration<3>::Instance()->SetType(LINEAR);
				WntConcentration<3>::Instance()->SetCellPopulation(cell_population);
				WntConcentration<3>::Instance()->SetCryptLength(crypt_length);


				// Create simulation object
				OffLatticeSimulationWithMonoclonalStoppingEvent simulator(cell_population);

				simulator.SetOutputDivisionLocations(true);
				simulator.SetDt(1.0/200.0);
				simulator.SetSamplingTimestepMultiple(200);
				simulator.SetEndTime(time_to_steady_state);

				// Add Volume Tracking Modifier
				MAKE_PTR(VolumeTrackingModifier<3>, p_modifier);
				simulator.AddSimulationModifier(p_modifier);

				//Create output directory
				std::stringstream out;
				out << "Run_" << sim_index << "/PercentMutant_" << percent_mutant;
				std::string output_directory = "CryptInvasionProbability/" +  out.str();
				simulator.SetOutputDirectory(output_directory);

				// Create a force law and pass it to the simulation
				MAKE_PTR(DifferentialAdhesionSpringForce<3>, p_force);
				p_force->SetMeinekeSpringStiffness(30.0); //normally 15.0 but 30 in all CellBased Papers;
				p_force->SetCutOffLength(cell_radius*3.0);
				simulator.AddForce(p_force);

				// Apply a retainer to keep stem and paneth cells at the base of the crypt
				MAKE_PTR(CellRetainerForce<3>, p_retainer_force);
				p_retainer_force->SetStemCellForceMagnitudeParameter(stem_retainer_force_magnitude);
				p_retainer_force->SetPanethCellForceMagnitudeParameter(paneth_retainer_force_magnitude);
				simulator.AddForce(p_retainer_force);

				// Apply a boundary condition to represent a 3d crypt
				MAKE_PTR_ARGS(CryptGeometryBoundaryCondition3d<3>, p_boundary_condition, (&cell_population, 0.0));
				simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

				// Create cell killer and pass in to crypt simulation
				MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer,(&cell_population, crypt_length*unit_vector<double>(3,2), unit_vector<double>(3,2)));
				simulator.AddCellKiller(p_cell_killer);

				// Run simulation
				simulator.Solve();



				// Now make a proportion of the cells mutant.

				// Iterate over all cells, to randomly assign a proportion to be mutant.
				for (AbstractCellPopulation<3>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				cell_iter != simulator.rGetCellPopulation().End();
				++cell_iter)
				{
				    // Generate a uniform random number to choose between Healthy and mutant cell in appropriate ratio

				    // Contact Inhibition specific parameters
				    if(!cell_iter->GetCellProliferativeType()->IsType<PanethCellProliferativeType>())
					{
						double u = RandomNumberGenerator::Instance()->ranf();
						if (u < percent_mutant)// Mutant cell
						{
							cell_iter->SetMutationState(p_mutant_state);
							// these cells are now mutant cells
							dynamic_cast<CryptCellCycleModel*>(cell_iter->GetCellCycleModel())->SetCellProliferationModel(mutant_cell_proliferation_model);
							dynamic_cast<CryptCellCycleModel*>(cell_iter->GetCellCycleModel())->SetIsWntDependentCellCycleDuration(mutant_wnt_dependend_ccd);

							dynamic_cast<CryptCellCycleModel*>(cell_iter->GetCellCycleModel())->SetWntThreshold(mutant_wnt_thresh);
							dynamic_cast<CryptCellCycleModel*>(cell_iter->GetCellCycleModel())->SetMutantWntThreshold(mutant_wnt_thresh);
							dynamic_cast<CryptCellCycleModel*>(cell_iter->GetCellCycleModel())->SetQuiescentVolumeFraction(mutant_CI);
						}
					}
				}

				simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();


				// Run the simulation to see the evolution of mutant v healthy cells.
				simulator.SetEndTime(time_to_steady_state + time_after_mutations);
				simulator.Solve();

				// Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
				WntConcentration<3>::Instance()->Destroy();
				SimulationTime::Instance()->Destroy();
				SimulationTime::Instance()->SetStartTime(0.0);
			}
			std::cout << " Runs complete \n" << std::flush;
		}
	}
};

#endif /*TESTCRYPTTAKEOVERPROBABILITYLITERATEPAPER_HPP_*/
