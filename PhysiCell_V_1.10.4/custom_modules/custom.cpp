/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
    
    Cell_Definition* pTumor_Up = find_cell_definition("Tumor_Up");
    Cell_Definition* pTumor_Down = find_cell_definition("Tumor_Down");
    Cell_Definition* pFibroblast = find_cell_definition("Fibroblast");
    Cell_Definition* pMacrophage = find_cell_definition("Macrophage");
    Cell_Definition* pEndothelial = find_cell_definition("Endothelial");
    Cell_Definition* pLung = find_cell_definition("Lung");
    
    pFibroblast->functions.update_phenotype = fibroblasts_phenotype_function;
    pEndothelial->functions.update_phenotype = endothelial_phenotype_function;
    pMacrophage->functions.update_phenotype = macrophages_phenotype_function; 
    pTumor_Down->functions.update_phenotype = tumor_down_phenotype_function; 
    pTumor_Up->functions.update_phenotype = tumor_up_phenotype_function; 
    pLung->functions.update_phenotype = phenotype_function;
    
     
    
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{	
	Cell* pC;
    Cell_Definition* pCD = find_cell_definition("Tumor_Down");
	pC = create_cell( *pCD ); 
    pC->assign_position( {525,510,0} );
    
	// load cells from your CSV file (if enabled)
	load_cells_csv(".\\config\\cells.csv"); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ 
    std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
    static Cell_Definition* pTumor_Up = find_cell_definition("Tumor_Up");
    static Cell_Definition* pTumor_Down = find_cell_definition("Tumor_Down");
    static Cell_Definition* pFibroblast = find_cell_definition("Fibroblast");
    static Cell_Definition* pMacrophage = find_cell_definition("Macrophage");
    static Cell_Definition* pEndothelial = find_cell_definition("Endothelial");
    static Cell_Definition* pLung = find_cell_definition("Lung");
    
		
    // Tumor UP -> black 
	if(pCell->type == pTumor_Up->type && pCell->custom_data["activated"]==0)
	{
		 output[0] = "rgb(0,0,0)"; 
		 output[2] = "rgb(0,0,0)"; 
	} else if (pCell->type == pTumor_Up->type && pCell->custom_data["activated"]==1){
        output[0] = "rgb(0,0,0)"; 
        output[2] = "rgb(0,0,0)"; 
    } 
    // Tumor down -> grey 
    else if (pCell->type == pTumor_Down->type && pCell->custom_data["activated"]==0){
        output[0] = "rgb(100,100,100)"; 
        output[2] = "rgb(100,100,100)"; 
    } else if (pCell->type == pTumor_Down->type && pCell->custom_data["activated"]==1){
        output[0] = "rgb(100,100,100)"; 
        output[2] = "rgb(100,100,100)"; 
    } 
    // Macrophage -> green 
    else if (pCell->type == pMacrophage->type && pCell->custom_data["activated"]==0){
        output[0] = "rgb(0,128,0)"; 
        output[2] = "rgb(0,128,0)"; 
    } else if (pCell->type == pMacrophage->type && pCell->custom_data["activated"]==1){
        output[0] = "rgb(0,110,0)"; 
        output[2] = "rgb(0,110,0)"; 
    } 
    // Endothelial -> tomato
    else if (pCell->type == pEndothelial->type && pCell->custom_data["activated"]==0){
        output[0] = "rgb(255,99,71)"; 
        output[2] = "rgb(255,99,71)"; 
    } else if (pCell->type == pEndothelial->type && pCell->custom_data["activated"]==1){
        output[0] = "rgb(235,79,51)"; 
        output[2] = "rgb(235,79,51)"; 
    } 
    // Lung -> violet 
    else if (pCell->type == pLung->type && pCell->custom_data["activated"]==0){
        output[0] = "rgb(238,130,238)"; 
        output[2] = "rgb(238,130,238)"; 
    } else if (pCell->type == pLung->type && pCell->custom_data["activated"]==1){
        output[0] = "rgb(238,130,238)"; 
        output[2] = "rgb(238,130,238)"; 
    } 
    // Fibroblast -> blueviolet, indigo
    else if (pCell->type == pFibroblast->type && pCell->custom_data["activated"]==0){
        output[0] = "rgb(138,43,226)"; 
        output[2] = "rgb(138,43,226)"; 
    } else if (pCell->type == pFibroblast->type && pCell->custom_data["activated"]==1){
        output[0] = "rgb(75,0,130)"; 
        output[2] = "rgb(75,0,130)"; 
    }
		
	return output;  }

	
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void activate_JNK_inhibition( void )
{
    static int JNK_inhibitor = microenvironment.find_density_index( "JNK_inhibitor" );
    microenvironment.set_substrate_dirichlet_activation(JNK_inhibitor,true);
}
void deactivate_JNK_inhibition( void )
{
    static int JNK_inhibitor = microenvironment.find_density_index( "JNK_inhibitor" );
    microenvironment.set_substrate_dirichlet_activation(JNK_inhibitor,true);
}
void activate_therapy( void )
{
    static int chemo = microenvironment.find_density_index( "chemo" );
    microenvironment.set_substrate_dirichlet_activation(chemo,true);
}
void deactivate_therapy( void )
{
    static int chemo = microenvironment.find_density_index( "chemo" );
    microenvironment.set_substrate_dirichlet_activation(chemo,true);
}


void custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void tumor_up_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{     
    
    // look for E_to_T -> increase growth rate or whatever 
    
    static int E_to_T_index = microenvironment.find_density_index( "E_to_T" );
    static int F_to_T_index = microenvironment.find_density_index( "F_to_T" );
    double chemokine_threshold = 1.0;
    double E_to_T_Concentration = (pCell->nearest_density_vector())[E_to_T_index];
    double F_to_T_Concentration = (pCell->nearest_density_vector())[F_to_T_index];
    static int apoptosis_model_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    
    if (E_to_T_Concentration > chemokine_threshold) {
            // Increase proliferation 
            pCell->phenotype.cycle.data.transition_rate(0,1) += 0.05*pCell->phenotype.cycle.data.transition_rate(0,1);
            // Decrease apoptosis 
            pCell->phenotype.death.rates[apoptosis_model_index] -= 0.05*pCell->phenotype.death.rates[apoptosis_model_index];
             
    }
    // look for F_to_T -> decrease apoptosis or whatever 
    if (F_to_T_Concentration > chemokine_threshold) {
            // Increase proliferation 
            pCell->phenotype.cycle.data.transition_rate(0,1) += 0.01*pCell->phenotype.cycle.data.transition_rate(0,1);
                         
    }
    
    // eat lung cells 
    // Get lung cell definitions
    static Cell_Definition* pLung = find_cell_definition("Lung");
    
    // Get neighborhood and see who is there 
    std::vector<Cell*> nearby = get_possible_neighbors( pCell); 
       
    for( int i=0 ; i < nearby.size() ; i++ )
    {
        Cell* pC = nearby[i]; 
        // if lung cell in the neighborhood, increase its apoptosis a lot
        std::vector<double> my_pos = pCell->position;
        std::vector<double> neighbor_pos = pC->position;
        double dist;
        dist = std::sqrt( (my_pos[0]-neighbor_pos[0])*(my_pos[0]-neighbor_pos[0]) + (my_pos[1]-neighbor_pos[1])*(my_pos[1]-neighbor_pos[1]) ); 
        if( pC->type == pLung->type && dist < 18.0) {
            if (pC->phenotype.death.rates[apoptosis_model_index] < 1e-8) {
                pC->phenotype.death.rates[apoptosis_model_index] = 1e-6;
            } else {
                pC->phenotype.death.rates[apoptosis_model_index] += 10*pC->phenotype.death.rates[apoptosis_model_index];
                break; 
            }
        } 
    }     
    
    
return; 
    
}

void tumor_down_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{     
            
    // eat lung cells 
    // Get lung cell definitions
    static Cell_Definition* pLung = find_cell_definition("Lung");
    static int apoptosis_model_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
    
    // Get neighborhood and see who is there 
    std::vector<Cell*> nearby = get_possible_neighbors( pCell); 
       
    for( int i=0 ; i < nearby.size() ; i++ )
    {
        Cell* pC = nearby[i]; 
        // if lung cell in the neighborhood, increase its apoptosis 100-fold
        std::vector<double> my_pos = pCell->position;
        std::vector<double> neighbor_pos = pC->position;
        double dist;
        dist = std::sqrt( (my_pos[0]-neighbor_pos[0])*(my_pos[0]-neighbor_pos[0]) + (my_pos[1]-neighbor_pos[1])*(my_pos[1]-neighbor_pos[1]) ); 
        if( pC->type == pLung->type && dist < 18.0) {
            if (pC->phenotype.death.rates[apoptosis_model_index] < 1e-8) {
                pC->phenotype.death.rates[apoptosis_model_index] = 1e-6;
            } else {
                pC->phenotype.death.rates[apoptosis_model_index] += 10*pC->phenotype.death.rates[apoptosis_model_index];
                break; 
            }
        }
    }        
    // see if chemo is on or off
    // get chemo concentration 
    static int chemo = microenvironment.find_density_index( "chemo" );
    double chemoConcentration = (pCell->nearest_density_vector())[chemo];
    if (chemoConcentration>1) {
            // induce cell apoptosis . highly increase apoptosis rate 
            if (pCell->phenotype.death.rates[apoptosis_model_index] < 1e-8) {
                pCell->phenotype.death.rates[apoptosis_model_index] = 1e-6;
            } else {
                pCell->phenotype.death.rates[apoptosis_model_index] += 10*pCell->phenotype.death.rates[apoptosis_model_index];
            }            
    }
    
    
    
return; 
    
}

void tumor_up_rule( Cell* pCell, Phenotype& phenotype, double dt )
{ 
    // Get cell definitions
    static Cell_Definition* pTumor_Up = find_cell_definition("Tumor_Up");
    static Cell_Definition* pTumor_Down = find_cell_definition("Tumor_Down");
    static Cell_Definition* pFibroblast = find_cell_definition("Fibroblast");
    static Cell_Definition* pMacrophage = find_cell_definition("Macrophage");
    static Cell_Definition* pEndothelial = find_cell_definition("Endothelial");
    static Cell_Definition* pLung = find_cell_definition("Lung");
    // Get neighborhood and see who is there 
    std::vector<Cell*> nearby = get_possible_neighbors( pCell); 
    int number_of_good_cells = 0; 
    int number_of_bad_cells = 0; 
    
    for( int i=0 ; i < nearby.size() ; i++ )
    {
        Cell* pC = nearby[i]; 
        // Is it a good cell ? 
        if( pC->type == pLung->type) {
            number_of_bad_cells++;
        } else {
            number_of_good_cells++;
        }
    }    
    // check the fraction of good to bad cells 
    //if (number_of_good_cells > number_of_bad_cells) {
    if (number_of_bad_cells > 0) {
        // generate random number 
        double random_val = 0;	
        int tumor_down_type_ID = pTumor_Down->type; 
        //SeedRandom();
        random_val = UniformRandom();
        // with 20% percent probability set transformation rate of this cell to tumor down to infinity 
        if (random_val < 0.8) {
                
                double transformation_rate = 9e9;
                set_single_behavior(pCell,"transform to cell type "+std::to_string(tumor_down_type_ID),transformation_rate); 
        }

    }
        
    // see if JNK inhibitor is on or off
    // get chemo concentration 
    static int JNK_inhibitor = microenvironment.find_density_index( "JNK_inhibitor" );
    double JNK_inhibitor_Concentration = (pCell->nearest_density_vector())[JNK_inhibitor];
    if (JNK_inhibitor_Concentration>1) {
            // transform JNK+ to JNK- 
        // generate random number 
        double random_val = 0;	
        int tumor_down_type_ID = pTumor_Down->type; 
        //SeedRandom();
        random_val = UniformRandom();
        // with 20% percent probability set transformation rate of this cell to tumor down to infinity 
        if (random_val < 0.8) {
                
                double transformation_rate = 9e9;
                set_single_behavior(pCell,"transform to cell type "+std::to_string(tumor_down_type_ID),transformation_rate); 
        }
                       
    }
    
return; }

void tumor_down_rule( Cell* pCell, Phenotype& phenotype, double dt )
{ 
    
    // Get cell definitions
    static Cell_Definition* pTumor_Up = find_cell_definition("Tumor_Up");
    static Cell_Definition* pTumor_Down = find_cell_definition("Tumor_Down");
    static Cell_Definition* pFibroblast = find_cell_definition("Fibroblast");
    static Cell_Definition* pMacrophage = find_cell_definition("Macrophage");
    static Cell_Definition* pEndothelial = find_cell_definition("Endothelial");
    static Cell_Definition* pLung = find_cell_definition("Lung");
    // Get neighborhood and see who is there 
    std::vector<Cell*> nearby = get_possible_neighbors( pCell); 
    int number_of_good_cells = 0; 
    int number_of_bad_cells = 0; 
    
    for( int i=0 ; i < nearby.size() ; i++ )
    {
        Cell* pC = nearby[i]; 
        // Is it a good cell ? 
        if( pC->type == pLung->type) {
            number_of_bad_cells++;
        } else {
            number_of_good_cells++;
        }
    }    
    // check the fraction of good to bad cells 
    //if (number_of_good_cells > number_of_bad_cells) {
    if (number_of_bad_cells == 0) {
        // generate random number 
        double random_val = 0;	
        int tumor_up_type_ID = pTumor_Up->type; 
        //SeedRandom();
        random_val = UniformRandom();
        // with 20% percent probability set transformation rate of this cell to tumor down to infinity 
        if (random_val < 0.8) {
                
                double transformation_rate = 9e9;
                set_single_behavior(pCell,"transform to cell type "+std::to_string(tumor_up_type_ID),transformation_rate); 
        }

    }
    
return; }

void fibroblasts_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{    
    static int chemokine_index = microenvironment.find_density_index( "T_to_F" );
    double chemokine_threshold = 1.0;
    double chemokineConcentration = (pCell->nearest_density_vector())[chemokine_index];
    if (chemokineConcentration > chemokine_threshold) {
            // Activate 
            pCell->custom_data["activated"] = 1; 
    }
    // If activate secrete 
        if (pCell->custom_data["activated"] == 1){
            // Set secretion 
            std::string substrate_name = "F_to_T"; 
            double value =  100; // secretion rate
            std::string behavior = substrate_name+" secretion";                        
            set_single_behavior(pCell,behavior,value);
                  
        }
    
return; }


void macrophages_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    
    static int chemokine_index = microenvironment.find_density_index( "T_to_M" );
    double chemokine_threshold = 1.0;
    double chemokineConcentration = (pCell->nearest_density_vector())[chemokine_index];
    if (chemokineConcentration > chemokine_threshold) {
            // Activate 
            pCell->custom_data["activated"] = 1; 
    }
    // If activate secrete 
        if (pCell->custom_data["activated"] == 1){
            // Set secretion 
            std::string substrate_name = "M_to_E"; 
            double value =  100; // secretion rate
            std::string behavior = substrate_name+" secretion";                        
            set_single_behavior(pCell,behavior,value);
                  
        }
    
return; }

void endothelial_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{    
    static int chemokine_index = microenvironment.find_density_index( "M_to_E" );
    double chemokine_threshold = 1.0;
    double chemokineConcentration = (pCell->nearest_density_vector())[chemokine_index];
    if (chemokineConcentration > chemokine_threshold) {
            // Activate 
            pCell->custom_data["activated"] = 1; 
    }
    // If activate secrete 
        if (pCell->custom_data["activated"] == 1){
            // Set secretion 
            std::string substrate_name = "E_to_T"; 
            double value =  100; // secretion rate
            std::string behavior = substrate_name+" secretion";                        
            set_single_behavior(pCell,behavior,value);
                  
        }
    // Deactivate endothelial ??? 
    
return; }

std::vector<Cell*> get_possible_neighbors( Cell* pCell)
{
	std::vector<Cell*> neighbors = {};
	
	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
		
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }
	
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end
		= pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();
	
	for( neighbor_voxel_index = pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin(); neighbor_voxel_index!= neighbor_voxel_index_end; ++neighbor_voxel_index)
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
		continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	return neighbors;
}	

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
