#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

#include "globals.h"
#include "geometry.h"
#include "MicroRoughnessHelper.h"

using namespace std;

material::material()                             
{
  theMicroRoughnessTable = 0;
  maxMicroRoughnessTable = 0;
  theMicroRoughnessTransTable = 0;
  maxMicroRoughnessTransTable = 0;

  theta_i_min = 0;
  theta_i_max = pi/2;
  Emin = 0;
  Emax = 1000*1e-9;
  no_theta_i = 90;
  noE = 1000;
  AngNoTheta = 180;
  AngNoPhi = 180;
  AngCut = 0.01;	
}

material::~material()
{
  /* if (theMicroRoughnessTable)      delete[] theMicroRoughnessTable;
  if (maxMicroRoughnessTable)      delete[] maxMicroRoughnessTable;
  if (theMicroRoughnessTransTable) delete[] theMicroRoughnessTransTable;
  if (maxMicroRoughnessTransTable) delete[] maxMicroRoughnessTransTable;
  cout << "THE MATERIAL DESTRUCTOR WAS CALLED!" << '\n'; */
}

double* material::GetMicroRoughnessTable ()
{
  return theMicroRoughnessTable;
}

double* material::GetMicroRoughnessTransTable ()
{
  return theMicroRoughnessTransTable;
}

void material::
               LoadMicroRoughnessTables(double* pMicroRoughnessTable,
                                        double* pmaxMicroRoughnessTable,
                                        double* pMicroRoughnessTransTable,
                                        double* pmaxMicroRoughnessTransTable)
{
  theMicroRoughnessTable      = pMicroRoughnessTable;
  maxMicroRoughnessTable      = pmaxMicroRoughnessTable;
  theMicroRoughnessTransTable = pMicroRoughnessTransTable;
  maxMicroRoughnessTransTable = pmaxMicroRoughnessTransTable;
}

void material::InitMicroRoughnessTables()
{
  

  cout << "thetadim: " << no_theta_i << " , Edim: " << noE << endl;

  // If both dimensions of the lookup-table are non-trivial:
  // delete old tables if existing and allocate memory for new tables

  if (no_theta_i*noE > 0) {
     if (theMicroRoughnessTable) delete[] theMicroRoughnessTable;
     theMicroRoughnessTable = new double[no_theta_i*noE];
     if (maxMicroRoughnessTable) delete[] maxMicroRoughnessTable;
     maxMicroRoughnessTable = new double[no_theta_i*noE];
     if (theMicroRoughnessTransTable) delete[] theMicroRoughnessTransTable;
     theMicroRoughnessTransTable = new double[no_theta_i*noE];
     if (maxMicroRoughnessTransTable) delete[] maxMicroRoughnessTransTable;
     maxMicroRoughnessTransTable = new double[no_theta_i*noE];
  }
}

void material::ComputeMicroRoughnessTables()
{
// Reads the parameters for the mr-probability computation from the
// corresponding material properties and stores it in the appropriate
// variables

  double b2 = RMSRough*(1.e-9)*RMSRough*(1.e-9);
  double w2 = CorrelLength*(1.e-9)*CorrelLength*(1.e-9);


  // The Fermi potential was saved in neV by P.F.

  double fermipot = FermiReal*(1.e-9);

  //cout << "Fermipot: " << fermipot/(1.e-9*eV) << "neV" << endl;

  double theta_i, E;

  // Calculates the increment in theta_i in the lookup-table
  theta_i_step = (theta_i_max-theta_i_min)/(no_theta_i-1);

  //cout << "theta_i_step: " << theta_i_step << endl;

  // Calculates the increment in energy in the lookup-table
  E_step = (Emax-Emin)/(noE-1);

  // Runs the lookup-table memory allocation
  InitMicroRoughnessTables();

  int counter = 0;

  cout << "Calculating MicroRoughness Tables...";

  // Writes the mr-lookup-tables to files for immediate control

  std::ofstream dateir("MRrefl.dat",std::ios::out);
  std::ofstream dateit("MRtrans.dat",std::ios::out);

  cout << theMicroRoughnessTable << endl;

  for (theta_i=theta_i_min; theta_i<=theta_i_max+1e-6; theta_i+=theta_i_step) {
      // Calculation for each cell in the lookup-table
      for (E=Emin; E<=Emax; E+=E_step) {
          *(theMicroRoughnessTable+counter) =
                      MicroRoughnessHelper::GetInstance() ->
                      IntIplus(E, fermipot, theta_i, AngNoTheta, AngNoPhi,
                               b2, w2, maxMicroRoughnessTable+counter, AngCut);

	  *(theMicroRoughnessTransTable+counter) =
                      MicroRoughnessHelper::GetInstance() -> 
                      IntIminus(E, fermipot, theta_i, AngNoTheta, AngNoPhi,
                                b2, w2, maxMicroRoughnessTransTable+counter,
                                AngCut);

          dateir << *(theMicroRoughnessTable+counter)      << endl;
          dateit << *(theMicroRoughnessTransTable+counter) << endl;

          counter++;

          //cout << counter << endl;
      }
  }

  dateir.close();
  dateit.close();

  // Writes the mr-lookup-tables to files for immediate control

  std::ofstream dateic("MRcheck.dat",std::ios::out);
  std::ofstream dateimr("MRmaxrefl.dat",std::ios::out);
  std::ofstream dateimt("MRmaxtrans.dat",std::ios::out);

  for (theta_i=theta_i_min; theta_i<=theta_i_max+1e-6; theta_i+=theta_i_step) {
      for (E=Emin; E<=Emax; E+=E_step) {

          // tests the GetXXProbability functions by writing the entries
          // of the lookup tables to files

          dateic  << GetMRIntProbability(theta_i,E)      << endl;
          dateimr << GetMRMaxProbability(theta_i,E)      << endl;
          dateimt << GetMRMaxTransProbability(theta_i,E) << endl;
      }
  }

  dateic.close();
  dateimr.close();
  dateimt.close();
}

double material::
                  GetMRIntProbability(double theta_i, double Energy)
{
  if (!theMicroRoughnessTable) {
     cout << "Dont have theMicroRoughnessTable" << endl;
     return 0.;
  }

  // if theta_i or energy outside the range for which the lookup table is
  // calculated, the probability is set to zero

  //cout << "theta_i: " << theta_i/degree << "degree"
  //       << "theta_i_min: " << theta_i_min/degree << "degree"
  //       << "theta_i_max: " << theta_i_max/degree << "degree"
  //       << "Emin: " << Emin/(1.e-9*eV) << "neV"
  //       << "Emax: " << Emax/(1.e-9*eV) << "neV" << endl;

  if (theta_i<theta_i_min || theta_i>theta_i_max || Energy<Emin || Energy>Emax)
     return 0.;

  // Determines the nearest cell in the lookup table which contains
  // the probability

  int theta_i_pos = int((theta_i-theta_i_min)/theta_i_step+0.5);
  int E_pos = int((Energy-Emin)/E_step+0.5);

  // lookup table is onoEensional (1 row), energy is in rows,
  // theta_i in columns

  //cout << "E_pos: " << E_pos << " theta_i_pos: " << theta_i_pos << endl;
  //cout << "Probability: " << *(theMicroRoughnessTable+E_pos+theta_i_pos*noE) << endl;

  return *(theMicroRoughnessTable+E_pos+theta_i_pos*noE);
}

double material::
                  GetMRIntTransProbability(double theta_i, double Energy)
{
  if (!theMicroRoughnessTransTable) return 0.;

  // if theta_i or energy outside the range for which the lookup table
  // is calculated, the probability is set to zero

  if (theta_i<theta_i_min || theta_i>theta_i_max || Energy<Emin || Energy>Emax)
     return 0.;

  // Determines the nearest cell in the lookup table which contains
  // the probability

  int theta_i_pos = int((theta_i-theta_i_min)/theta_i_step+0.5);
  int E_pos = int((Energy-Emin)/E_step+0.5);

  // lookup table is onoEensional (1 row), energy is in rows,
  // theta_i in columns

  return *(theMicroRoughnessTransTable+E_pos+theta_i_pos*noE);
}

double material::
                  GetMRMaxProbability(double theta_i, double Energy)
{
  if (!maxMicroRoughnessTable) return 0.;

  // if theta_i or energy outside the range for which the lookup table
  // is calculated, the probability is set to zero

  if (theta_i<theta_i_min || theta_i>theta_i_max || Energy<Emin || Energy>Emax)
     return 0.;

  // Determines the nearest cell in the lookup table which contains
  // the probability

  int theta_i_pos = int((theta_i-theta_i_min)/theta_i_step+0.5);
  int E_pos = int((Energy-Emin)/E_step+0.5);

  // lookup table is onoEensional (1 row), energy is in rows,
  // theta_i in columns

  return *(maxMicroRoughnessTable+E_pos+theta_i_pos*noE);
}

double material::
                  GetMRMaxTransProbability(double theta_i, double Energy)
{
  if (!maxMicroRoughnessTransTable) return 0.;

  // if theta_i or energy outside the range for which the lookup table
  // is calculated, the probability is set to zero

  if (theta_i<theta_i_min || theta_i>theta_i_max || Energy<Emin || Energy>Emax)
     return 0.;

  // Determines the nearest cell in the lookup table which contains
  // the probability

  int theta_i_pos = int((theta_i-theta_i_min)/theta_i_step+0.5);
  int E_pos = int((Energy-Emin)/E_step+0.5);

  // lookup table is onoEensional (1 row), energy is in rows,
  // theta_i in columns

  return *(maxMicroRoughnessTransTable+E_pos+theta_i_pos*noE);
}

double material::
                  GetMRProbability(double theta_i, double Energy,
                                   double fermipot,
                                   double theta_o, double phi_o)
{
  return MicroRoughnessHelper::GetInstance()->
          ProbIplus(Energy, fermipot, theta_i, theta_o, phi_o, RMSRough, CorrelLength, AngCut);
}

double material::
                  GetMRTransProbability(double theta_i, double Energy,
                                        double fermipot,
                                        double theta_o, double phi_o)
{
  return MicroRoughnessHelper::GetInstance()->
          ProbIminus(Energy, fermipot,theta_i, theta_o, phi_o, RMSRough, CorrelLength, AngCut);
}

bool material::ConditionsValid(double E,
                                                     double VFermi,
                                                     double theta_i)
{
  double k =   sqrt(2*neutron_mass_c2*E      / hbarc_squared);
  double k_l = sqrt(2*neutron_mass_c2*VFermi / hbarc_squared);

  //cout << " Energy: " << E/(1.e-9*eV) << "neV"
  //       << " VFermi: " << VFermi/(1.e-9*eV) << "neV"
  //       << " theta:  " << theta_i/degree << "degree" << endl;

  //cout << " Conditions - 2*b*k*cos(theta): " << 2*b*k*cos(theta_i)
  //       << ", 2*b*kl: "                       << 2*b*k_l << endl;

  // see eq. 17 of the Steyerl paper

  if (2*RMSRough*k*cos(theta_i) < 1 && 2*RMSRough*k_l < 1) return true;
  else return false;
}

bool material::TransConditionsValid(double E,
                                                          double VFermi,
                                                          double theta_i)
{
  double k2   = 2*neutron_mass_c2*E      / hbarc_squared;
  double k_l2 = 2*neutron_mass_c2*VFermi / hbarc_squared;

  if (E*cos(theta_i)*cos(theta_i) < VFermi) return false;

  double kS2 = k_l2 - k2;

  //cout << "Conditions; 2bk' cos(theta): " << 2*b*sqrt(kS2)*cos(theta_i) << 
  //          ", 2bk_l: " << 2*b*sqrt(k_l2) << endl;

  // see eq. 18 of the Steyerl paper

  if (2*RMSRough*sqrt(kS2)*cos(theta_i) < 1 && 2*RMSRough*sqrt(k_l2) < 1) return true;
  else return false;
}


TGeometry::TGeometry(TConfig &geometryin){
	vector<material> materials;

	for (map<string, string>::iterator i = geometryin["MATERIALS"].begin(); i != geometryin["MATERIALS"].end(); i++){
		material mat;
		mat.name = i->first;
		istringstream ss(i->second);

		ss >> mat.FermiReal >> mat.FermiImag >> mat.DiffProb >> mat.SpinflipProb >> mat.RMSRough >> mat.CorrelLength >> mat.Interaction;
		
		if (mat.Interaction == 2) {
			cout << "ABOUT TO CALCULATE MR TABLES!!! for: " << mat.name << '\n';
			mat.ComputeMicroRoughnessTables();
			cout << "CALCULATED MR TABLES!!! WOOOOOOOOOOOOOOOOOOOO!" << '\n';
		}

		if (ss)
			materials.push_back(mat);
		else
			cout << "Could not load material " << i->first << '\n';

	
	}

	string line;
	string STLfile;
	string matname;
	char name[80];
	for (map<string, string>::iterator i = geometryin["GEOMETRY"].begin(); i != geometryin["GEOMETRY"].end(); i++){	// parse STLfile list
		solid model;
		istringstream(i->first) >> model.ID;
		if (model.ID < 0){
			cout << "You defined a solid with ID " << model.ID << " < 0! IDs have to be larger than zero!\n";
			exit(-1);
		}
		for (vector<solid>::iterator j = solids.begin(); j != solids.end(); j++){
			if (j->ID == model.ID){
				cout << "You defined solids with identical ID " << model.ID << "! IDs have to be unique!\n";
				exit(-1);
			}
		}

		istringstream ss(i->second);
		ss >> STLfile >> matname;
		if (ss){
			for (unsigned i = 0; i < materials.size(); i++){
				if (matname == materials[i].name){
					model.mat = materials[i];
					if (model.ID > 1){
						mesh.ReadFile(STLfile.c_str(),solids.size(),name);
						model.name = name;
					}
					else
						model.name = "default solid";
					break;
				}
				else if (i+1 == materials.size()){
					printf("Material %s used for %s but not defined in geometry.in!",matname.c_str(),name);
					exit(-1);
				}
			}
		}
		else{
			cout << "Could not load solid with ID " << model.ID << "! Did you define invalid parameters?\n";
			exit(-1);
		}

		long double ignorestart, ignoreend;
		while (ss){
			ss >> ignorestart;
			if (!ss) // no more ignore times found
				break;
			char dash;
			ss >> dash;
			ss >> ignoreend;
			if (ss && dash == '-'){
				model.ignoretimes.push_back(ignorestart);
				model.ignoretimes.push_back(ignoreend);
			}
			else{
				cout << "Invalid ignoretimes in geometry.in" << '\n';
				exit(-1);
			}
		}

		if (model.ID == 1)
			defaultsolid = model;
		else
			solids.push_back(model);
	}
	mesh.Init();
	cout << '\n';
}


bool TGeometry::GetCollisions(const long double x1, const long double p1[3], const long double h, const long double p2[3], map<TCollision, bool> &colls){
	set<TCollision> c;
	mesh.Collision(p1,p2,c);
	colls.clear();
	for (set<TCollision>::iterator it = c.begin(); it != c.end(); it++){
		colls[*it] = false;
		std::vector<long double> *times = &solids[it->sldindex].ignoretimes;
		if (!times->empty()){
			long double x = x1 + h*it->s;
			for (unsigned int i = 0; i < times->size(); i += 2){
				if (x >= (*times)[i] && x < (*times)[i+1]){
					colls[*it] = true; // set ignored flag if collision should be ignored according to geometry.in
					break;
				}
			}
		}
	}
	return !colls.empty();
}


void TGeometry::GetSolids(const long double t, const long double p[3], std::map<solid, bool> &currentsolids){
	long double p2[3] = {p[0], p[1], mesh.tree.bbox().zmin() - REFLECT_TOLERANCE};
	map<TCollision, bool> c;
	currentsolids.clear();
	currentsolids[defaultsolid] = false;
	if (GetCollisions(t,p,0,p2,c)){	// check for collisions of a vertical segment from p to lower border of bounding box
		for (map<TCollision, bool>::iterator i = c.begin(); i != c.end(); i++){
			solid *sld = &solids[i->first.sldindex];
			if (currentsolids.count(*sld) > 0) // if there is a collision with a solid already in the list, remove it from list
				currentsolids.erase(*sld);
			else
				currentsolids[*sld] = i->second; // else add solid to list
		}
	}
}

