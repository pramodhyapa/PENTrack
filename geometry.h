/**
 * \file
 * Contains class to include experiment geometry.
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <limits>
#include <map>

#include "trianglemesh.h"
#include "globals.h"


using namespace std;

static const double REFLECT_TOLERANCE = 1e-8;  ///< max distance of reflection point to actual surface collision point

/// Struct to store material properties (read from geometry.in, right now only for neutrons)
struct material{
	std::string name; ///< Material name
	long double FermiReal; ///< Real part of Fermi potential
	long double FermiImag; ///< Imaginary part of Fermi potential
	long double DiffProb; ///< Diffuse reflection probability for Lambert reflection
	long double SpinflipProb; ///< Probability for spin flip on reflection
	int Interaction;

	long double RMSRough; //<The root mean square roughness, b, from the MicroRoughness model
	long double CorrelLength; //<Correlation length, w, from the MicroRoughness model
	
	// Pointer to the integral reflection probability table
  	double* theMicroRoughnessTable;
  	// Pointer to the maximum reflection probability table
  	double* maxMicroRoughnessTable;
  	// Pointer to the integral transmission probability table
  	double* theMicroRoughnessTransTable;
  	// Pointer to the maximum transmission probability table
  	double* maxMicroRoughnessTransTable;

  	double theta_i_min;
  	double theta_i_max;
  	double Emin;
  	double Emax;
  	int no_theta_i;
  	int noE;
  	double theta_i_step;
  	double E_step;
	int AngNoTheta;
	int AngNoPhi;
  	double AngCut;
	
	material();
  	~material();

  	// returns the pointer to the mr-reflection table
  	double* GetMicroRoughnessTable();

  	// returns the pointer to the mr-transmission table
  	double* GetMicroRoughnessTransTable();

  	// Assigns double-array to the table-pointers, currently not used
  	void LoadMicroRoughnessTables(double*, double*, double*, double*);

  	// Creates new double arrays and assigns them to the table pointers
  	void InitMicroRoughnessTables();

  	// Reads the MR-parameters from the corresponding fields and starts
  	// the computation of the mr-tables
  	void ComputeMicroRoughnessTables();

  	// returns the integral prob. value for a theta_i - E pair
  	double GetMRIntProbability (double, double);

  	// returns the maximum prob. value for a theta_i - E pair
  	double GetMRMaxProbability (double, double);

  	// returns the mr-prob.

  	// arguments:
  	//         1) theta_i
  	//         2) Energy
  	//         3) V_F
  	//         4) theta_o
  	//         5) phi_o

  	double GetMRProbability (double, double, double, double, double);

  	// returns the integral transmission prob. value for a theta_i - E pair
  	double GetMRIntTransProbability (double, double);

  	// returns the maximum transmission prob. for a theta_i - E pair
  	double GetMRMaxTransProbability (double, double);

  	// returns the mr-transmission-prob.

  	// arguments:
  	//         1) theta_i
  	//         2) E
  	//         3) V_F
  	//         4) theta_o
  	//         5) phi_o

  	double GetMRTransProbability (double, double,
                                  double, double, double);

  	// Checks if the validity condition for the microroughness model are
  	// satisfied, cf. Steyerl-paper p. 175
  	bool ConditionsValid (double E, double VFermi, double theta_i);

  	// Checks if the validity conditions for the transmission of the
  	// microroughness model are satisfied
  	bool TransConditionsValid (double E, double VFermi, double theta_i);

};


/// Struct to store solid information (read from geometry.in)
struct solid{
	std::string name; ///< name of solid
	material mat; ///< material of solid
	unsigned ID; ///< ID of solid
	std::vector<long double> ignoretimes; ///< pairs of times, between which the solid should be ignored

	/**
	 * Comparison operator used to sort solids by priority (descending)
	 *
	 * @param s solid struct compared to this solid
	 *
	 * @return Returns true if this solid's ID is larger than the other one's
	 */
	bool operator< (const solid s) const { return ID > s.ID; };
};


/**
 * Class to include experiment geometry.
 *
 * Loads solids and materials from geometry.in, maintains solids list, checks for collisions.
 */
struct TGeometry{
	public:
		TTriangleMesh mesh; ///< kd-tree structure containing triangle meshes from STL-files
		vector<solid> solids; ///< solids list
		solid defaultsolid; ///< "vacuum", this solid's properties are used when the particle is not inside any other solid
		
		/**
		 * Constructor, reads geometry configuration file, loads triangle meshes.
		 *
		 * @param geometryin TConfig struct containing MATERIALS and GEOMETRY config section
		 */
		TGeometry(TConfig &geometryin);


		/**
		 * Check if point is inside geometry bounding box.
		 *
		 * @param y1 Position vector of segment start
		 * @param y2 Position vector of segment end
		 *
		 * @return Returns true if point is inside the bounding box
		 */
		bool CheckSegment(const long double y1[3], const long double y2[3]){
			return CGAL::do_intersect(mesh.tree.bbox(), CSegment(CPoint(y1[0], y1[1], y1[2]), CPoint(y2[0], y2[1], y2[2])));
		};
		

		/**
		 * Checks if line segment p1->p2 collides with a surface.
		 *
		 * Calls KDTree::Collision to check for collisions and removes all collisions
		 * which should be ignored (given by ignore times in geometry configuration file).
		 *
		 * @param x1 Start time of line segment
		 * @param p1 Start point of line segment
		 * @param h Time length of line segment
		 * @param p2 End point of line segment
		 * @param colls List of collisions, paired with bool indicator it it should be ignored
		 *
		 * @return Returns true if line segment collides with a surface
		 */
		bool GetCollisions(const long double x1, const long double p1[3], const long double h, const long double p2[3], map<TCollision, bool> &colls);
		
			
		/**
		 * Get solids in which the point p lies
		 *
		 * @param t Time
		 * @param p Point to test
		 * @param currentsolids Map of solids in which the point is inside paired with information if it was ignored or not
		 */
		void GetSolids(const long double t, const long double p[3], map<solid, bool> &currentsolids);
};

#endif /*GEOMETRY_H_*/
