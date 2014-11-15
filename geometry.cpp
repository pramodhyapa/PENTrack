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


using namespace std;

material::material()                             
{
  AngNoTheta = 180;
  AngNoPhi = 180;
  AngCut = 0.01;	
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

