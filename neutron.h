/**
 * \file
 * Neutron class definition.
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

static const char* NAME_NEUTRON = "neutron";

#include <complex>

#include "globals.h"
#include "particle.h"
#include "source.h"
#include "proton.h"
#include "electron.h"
#include "MRHelperPointerless.h"

static const int MAX_DICE_ROLL = 42000000; ///< number of tries to find neutron start point

/**
 * Neutron particle class.
 *
 * Simulates an ultra cold neutron including gravitation, magnetic forces on its magnetic dipole moment
 * and tries to estimate spin flip probability in magnetic fields.
 */
struct TNeutron: TParticle{
public:
	/**
	 * Constructor, create neutron, set start values randomly according to particle.in.
	 *
	 * Sets start time TParticle::tstart according to [SOURCE] in geometry.in.
	 * Sets start TParticle::polarisation according to particle.in.
	 * Sets start energy according to TMCGenerator::NeutronSpectrum.
	 * Then tries to find a position according to [SOURCE] in geometry.in.
	 * For volume sources the total energy TParticle::Hstart is diced by TMCGenerator::NeutronSpectrum
	 * and the position search is repeated until the kinetic energy TParticle::Estart is >0.
	 * Additionally the kinetic energy spectrum is weighted by sqrt(Estart) (see Golub/Richardson/Lamoreaux p. 82).
	 * For surface sources the kinetic energy TParticle::Estart is diced and all positions are allowed.
	 *
	 * @param number Particle number
	 * @param ageometry Geomtry in which the particle will be simulated
	 * @param src TSource in which particle should be generated
	 * @param mcgen TMCGenerator used to dice inital values
	 * @param afield TFieldManager used to calculate energies
	 */
	TNeutron(int number, TGeometry &ageometry, TSource &src, TMCGenerator &mcgen, TFieldManager *afield)
			: TParticle(NAME_NEUTRON, 0, m_n, mu_nSI, gamma_n){
		long double t;
		int polarisation = mcgen.DicePolarisation(name);
		long double H = mcgen.Spectrum(name);
		long double E;
		long double p[3];

		cout << "Dice starting position for E_neutron = " << H*1e9 << " neV ";
		long double phi, theta;
		for (int nroll = 0; nroll <= MAX_DICE_ROLL; nroll++){ // try to create particle only MAX_DICE_ROLL times
			if (nroll % 1000000 == 0){
				printf("."); // print progress
			}
			src.RandomPointInSourceVolume(mcgen, t, H, p[0], p[1], p[2], phi, theta);
			ageometry.GetSolids(t, p, currentsolids);
			if (src.sourcemode == "volume" || src.sourcemode == "customvol"){ // dice H for volume source
				E = H - Epot(t, p, polarisation, afield, currentsolids);
				if (mcgen.UniformDist(0,1) > sqrt(E/H))
					E = -1; // correct weighting of phase space (see Golub/Richardson/Lamoreaux p. 82)
			}
			else{ // dice E for surface source
				E = H;
				H = E + Epot(t, p, polarisation, afield, currentsolids);
			}
			if (E >= 0){
				printf(" Found! %i dice(s) rolled\n", nroll+1);
				break;
			}
			else if (nroll >= MAX_DICE_ROLL){
				ID = ID_INITIAL_NOT_FOUND;
				printf(" ABORT: Failed %i times to find a compatible spot!! NO particle will be simulated!!\n\n", MAX_DICE_ROLL);
				return;
			}
		}
		InitE(number, t, mcgen.LifeTime(name), p[0], p[1], p[2],
			E, phi, theta, polarisation, mcgen.MaxTrajLength(name), ageometry, afield);
	};



protected:
	/**
	 * Check for reflection on surfaces.
	 *
	 * Uses Fermi-potential formalism to calculate reflection/transmission probabilities and uses MicroRoughness enhancement if that model is used.
	 * Depending on the Interaction model, each reflection has a certain chance of being diffuse/MicroRoughness diffuse.
	 * Diffuse reflection angles are cosine-distributed around the normal vector of the surface.
	 * Logs surface hits if reflectlog is set in config.in.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param polarisation Polarisation of particle, may be altered
	 * @param normal Normal vector of hit surface
	 * @param leaving Solid that the particle is leaving
	 * @param entering Solid that the particle is entering (can be modified by method)
	 * @param trajectoryaltered Returns true if the particle trajectory was altered
	 * @param traversed Returns true if the material boundary was traversed by the particle
	 */
	void OnHit(long double x1, long double y1[6], long double &x2, long double y2[6], int &polarisation,
				const long double normal[3], solid *leaving, solid *entering, bool &trajectoryaltered, bool &traversed){
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
		long double E = Ekin(&y1[3]); //Kinetic Energy of Neutron REDUNDANT POSSIBLY
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
		long double VFentering = entering->mat.FermiReal*1e-9; //Real component of the Fermi potential of entering material in [eV]
		long double VFleaving = leaving->mat.FermiReal*1e-9; //Real component of the Fermi potential of leaving material in [eV]
		int Interaction = entering->mat.Interaction; //This determines whether it is specular, difuse or MicroRoughness reflection/transmission
		MRHelperPointerless MRFunctions;

		//cout << "This is the energy of the neutron: " << E << endl;

		bool refl = false;
		trajectoryaltered = false;
		traversed = true;

		long double prob = mc->UniformDist(0,1);
		long double theta_i = IncidentTheta(y1, normal); //the incident neutron angle in radians
		long double Estep = VFentering - VFleaving;

		//cout << "This is the Incident Angle: " << theta_i/conv << '\n' ;

		//========= TRANSMISSION ========//
		if (Enormal > Estep){ // transmission only possible if E > Estep
			//*******************
			long double transprob = TransmissionProb(Enormal, Estep, entering, Interaction);
			if (prob < transprob){ // -> transmission
				switch (Interaction) { //This switch changes to the appropriate model depending on what the int Interaction is
 					case 0: { //Specular case
						//************** Snell's Law Transmission **************
						DoSnellTransmission(y1, y2, normal, Estep);
					}
					break;
  					case 1: { //Specular-diffuse case
						if (prob >= entering->mat.DiffProb){
							DoSnellTransmission(y1, y2, normal, Estep);
						}
						else{
							//************** Diffuse transmission ************
							DoDiffuseTransmission(x1, y1, x2, y2, normal, Estep);
						}
					}
					break;
					case 2: { //Specular-MicroRoughness case
						bool valid = entering->mat.TransConditionsValid(E, Estep, theta_i);
						if (valid) { //Do MR transmission only if the conditions are valid
							cout << "MR Transmission conditions are valid!" << endl;
							int AngNoTheta = entering->mat.AngNoTheta;
							int AngNoPhi = entering->mat.AngNoPhi;
							double AngCut = entering->mat.AngCut;
							long double b = entering->mat.RMSRough*1e-9; // RMS roughness of entering material in [m]
							long double w = entering->mat.CorrelLength*1e-9; // Correlation length of entering material in [m]
							long double prob1 = mc->UniformDist(0,1); //random dice to determine if it is MRM or Snell's Law
							long double MRprob = MRFunctions.IntIminus(E, Estep, theta_i, AngNoTheta, AngNoPhi, b*b, w*w, AngCut);
							//cout << '\n' <<'\n' << " MicroRoughness Prob = " << diffuseprob << '\n' << '\n' << '\n';

							if (prob1<MRprob) { //this is for MRM transmission
								//************** MicroRoughness Transmission **************
								DoMicroRoughnessTransmission(Estep, x1, y1, x2, y2, normal, entering, MRFunctions);
								traversed = true;
								trajectoryaltered = true;
							}

							else {
								//************** Snell's Law Transmission **************
								DoSnellTransmission(y1, y2, normal, Estep);
							}

						}
						else { //If the MR transmission condition is not valid, switch to Diffuse-specular
							if (prob >= entering->mat.DiffProb){
								//************** Snell's Law Transmission **************
								DoSnellTransmission(y1, y2, normal, Estep);
							}
							else{
								//************** Diffuse transmission ************
								DoDiffuseTransmission(x1, y1, x2, y2, normal, Estep);
							}
						}
					}
     					break;
  					default:
    					cout << "Error in determining type of interaction for transmission!" << endl;
					break;
				}
			}
			else // no transmission -> reflection
				refl = true;
			//*******************

		}

		//======== ABSORPTION =====///
		else{
			long double reflprob = ReflectionProb(Enormal, Estep, entering, Interaction);
			if (prob > reflprob){ // -> absorption on reflection
				DoAbsorption(x1, y1, x2, y2, entering);
				traversed = false;
				trajectoryaltered = true;
			}
			else // no absorption -> reflection
				refl = true;
		}

		//=========== REFLECTION =========//
		if (refl){
			switch (Interaction) { //This switch changes to the appropriate model depending on what value the int Interaction takes
				case 0: { //Specular case
					//************** Specular Reflection **************
					DoSpecularReflection(normal, x1, y1, x2, y2);
				}
				break;
 					case 1: { //Specular-diffuse case
					if (prob >= entering->mat.DiffProb){
						//************** Specular Reflection **************
						DoSpecularReflection(normal, x1, y1, x2, y2);
					}
					else{
						//************** Diffuse reflection ************
						DoDiffuseReflection(x1, y1, x2, y2, normal);
					}
				}
				break;
				case 2: { //Specular-MicroRoughness case
					bool valid = entering->mat.ConditionsValid(E, Estep, theta_i);
					if (valid) { //Do MR transmission only if the conditions are valid
					//cout << "MR Transmission conditions are valid!" << endl;
						int AngNoTheta = entering->mat.AngNoTheta;
						int AngNoPhi = entering->mat.AngNoPhi;
						double AngCut = entering->mat.AngCut;
						long double b = entering->mat.RMSRough*1e-9; // RMS roughness of entering material in [m]
						long double w = entering->mat.CorrelLength*1e-9; // Correlation length of entering material in [m]
						long double prob1 = mc->UniformDist(0,1); //random dice to determine if it is MRM or Snell's Law
						long double MRprob = MRFunctions.IntIplus(E, Estep, theta_i, AngNoTheta, AngNoPhi, b*b, w*w, AngCut);
						//cout << '\n' <<'\n' << " MicroRoughness Prob = " << diffuseprob << '\n' << '\n' << '\n';

						if (prob1<MRprob) { //this is for MRM reflection
							//************** MicroRoughness Reflection **************
							DoMicroRoughnessReflection(Estep, x1, y1, x2, y2, normal, entering, MRFunctions);
						}

						else {
							//************** Specular Reflection **************
							DoSpecularReflection(normal, x1, y1, x2, y2);
						}

					}
					else { //If the MR reflection condition is not valid, switch to Diffuse-specular
						if (prob >= entering->mat.DiffProb){
							//************** Specular Reflection **************
							DoSpecularReflection(normal, x1, y1, x2, y2);
						}
						else{
							//************** Diffuse reflection ************
							DoDiffuseReflection(x1, y1, x2, y2, normal);
						}
					}
				}
				break;
				default:
				cout << "Error in determining type of interaction for reflection!" << endl;
				break;
			}

			if (mc->UniformDist(0,1) < entering->mat.SpinflipProb){
				polarisation *= -1;
				Nspinflip++;
			}
			
			trajectoryaltered = true;
			traversed = false;
		}
	};


	/**
	 * Checks for absorption in solids using Fermi-potential formalism and does some additional calculations for neutrons
	 *
	 * Caclulates adiabacity condition (TNeutron::frac).
	 * Estimates spin flip probability according to Vladimirskii (TNeutron::vlad)
	 * and by doing a bruteforce integration of the Bloch equation describing spin precession.
	 * Additionally it can print out a spatial neutron distribution matrix.
	 *
	 * @param x1 Start time of line segment
	 * @param y1 Start point of line segment
	 * @param x2 End time of line segment, may be altered
	 * @param y2 End point of line segment, may be altered
	 * @param currentsolid Solid through which the particle is moving
	 * @return Returns true if particle was absorbed
	 */
	bool OnStep(long double x1, long double y1[6], long double &x2, long double y2[6], solid currentsolid){
		bool result = false;
		if (currentsolid.mat.FermiImag > 0){
			long double prob = mc->UniformDist(0,1);
			complex<long double> E(0.5*m_n*(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]), currentsolid.mat.FermiImag*1e-9); // E + i*W
			complex<long double> k = sqrt(2*m_n*E)*ele_e/hbar; // wave vector
			long double l = sqrt(pow(y2[0] - y1[0], 2) + pow(y2[1] - y1[1], 2) + pow(y2[2] - y1[2], 2)); // travelled length
			if (prob > exp(-imag(k)*l)){ // exponential probability decay
				x2 = x1 + mc->UniformDist(0,1)*(x2 - x1); // if absorbed, chose a random time between x1 and x2
				for (int i = 0; i < 6; i++)
					y2[i] = stepper->dense_out(i, x2, stepper->hdid); // set y2 to position at random time
				ID = currentsolid.ID;
				printf("Absorption!\n");
				result = true; // stop integration
			}
		}

		// do special calculations for neutrons (spinflipcheck, snapshots, etc)
		if (neutdist == 1)
			fillndist(x1, y1, x2, y2); // write spatial neutron distribution
/*
		if (field){
			long double B[4][4];
			field->BField(y1[0],y1[1],y1[2],x1,B);
*//*
			if (B[3][0] > 0){
				// spin flip properties according to Vladimirsky and thumbrule
				vlad = vladimirsky(B[0][0], B[1][0], B[2][0],
								   B[0][1], B[0][2], B[0][3], B[1][1], B[1][2], B[1][3], B[2][1], B[2][2], B[2][3], B[3][0],
								   y1[3], y1[4], y1[5]);
				frac = thumbrule(B[0][0], B[1][0], B[2][0],
								   B[0][1], B[0][2], B[0][3], B[1][1], B[1][2], B[1][3], B[2][1], B[2][2], B[2][3], B[3][0],
								   y1[3], y1[4], y1[5]);
				vladtotal *= 1-vlad;
				if (vlad > 1e-99)
					vladmax = max(vladmax,log10(vlad));
				if (frac > 1e-99)
					thumbmax = max(thumbmax,log10(frac));
			}
*/
/*			long double B2[4][4];
			field->BField(y2[0],y2[1],y2[2],x2,B2);
			long double sp = BruteForceIntegration(x1,y1,B,x2,y2,B2); // integrate spinflip probability
//			if (1-sp > 1e-30) logBF = log10(1-sp);
//			else logBF = -99;
//			BFsurvprob *= sp;
			// flip the spin with a probability of 1-BFsurvprob
			if (flipspin && mc->UniformDist(0,1) < 1-sp)
			{
				polarisation *= -1;
				Nspinflip++;
				printf("\n The spin has flipped! Number of flips: %i\n",Nspinflip);
				result = true;
			}
		}
*/
		return result;
	};


	/**
	 * Neutron decay.
	 *
	 * Create decay proton and electron and add them to secondaries list
	 */
	void Decay(){
		long double v_p[3], v_e[3];
		TParticle *p;
		mc->NeutronDecay(&yend[3], v_p, v_e);
		p = new TProton(particlenumber, tend, mc->LifeTime("proton"),
							yend[0], yend[1], yend[2], v_p[0], v_p[1], v_p[2],
							mc->DicePolarisation("proton"), mc->MaxTrajLength("proton"), *geom, field);
		secondaries.push_back(p);
		p = new TElectron(particlenumber, tend, mc->LifeTime("electron"),
							yend[0], yend[1], yend[2], v_e[0], v_e[1], v_e[2],
							mc->DicePolarisation("electron"), mc->MaxTrajLength("electron"), *geom, field);
		secondaries.push_back(p);
	};


	/**
	 * Potential energy, other than electromagnetic or gravitational, which is added to total energy of neutron
	 *
	 * @param t Time
	 * @param y Position vector
	 * @param polarisation Particle polarisation
	 * @param field TFieldManager for electromagnetic potential
	 * @param solids List of solids for optional material potential
	 *
	 * @return Returns potential energy plus Fermi-Potential of solid
	 */
	long double Epot(const long double t, const long double y[3], int polarisation, TFieldManager *field, map<solid, bool> &solids){
		map<solid, bool>::iterator sld = solids.begin();
		while (sld->second)
			sld++;
		return TParticle::Epot(t, y, polarisation, field, solids) + sld->first.mat.FermiReal*1e-9;
	}

	long double TransmissionProb(long double Enormal, long double Estep, solid *entering, int I) { //This returns the probability of transmission, modified by the MR model if the int I == 2
		long double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
		long double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
		long double transprob = 4*k1*k2/(k1 + k2)/(k1 + k2); // transmission probability
		cout << " This is the transmission probability before MR modification = " << transprob << '\n';
		if (I == 2) {
			long double b = entering->mat.RMSRough*1e-9; // RMS roughness of entering material in [m]
			long double w = entering->mat.CorrelLength*1e-9; // Correlation length of entering material in [m]
			long double Vr = entering->mat.FermiReal*1e-9; //Real component of the Fermi potential of entering material in [eV]
			long double kl = sqrt(2*neutron_mass_c2*Vr)/hbarc_squared; //limiting wave number defined in Heule's thesis after Eq 4.11
			transprob = transprob*sqrt(1+(2*b*b*kl*kl)/(1+0.85*kl*w+2*kl*kl*w*w)); // transmission probability from Eq 4.26 of Heule's thesis
			cout << " This is the transmission probability after MR modification = " << transprob << '\n';
		}
		return transprob;

	}

	long double ReflectionProb(long double Enormal, long double Estep, solid *entering, int I) { //This returns the probability of transmission, modified by the MR model if the int I == 2
		long double k1 = sqrt(Enormal); // wavenumber in first solid (only real part)
		complex<long double> iEstep(Estep, -entering->mat.FermiImag*1e-9); // potential step using imaginary potential V - i*W of second solid
		complex<long double> k2 = sqrt(Enormal - iEstep); // wavenumber in second solid (including imaginary part)
		long double reflprob = pow(abs((k1 - k2)/(k1 + k2)), 2); // reflection probability
		//cout << " This is the reflection probability before MR modification = " << reflprob << '\n';
		if (I == 2) {
			long double k3 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for transmission!)
			long double b = entering->mat.RMSRough*1e-9; // RMS roughness of entering material in [m]
			long double w = entering->mat.CorrelLength*1e-9; // Correlation length of entering material in [m]
			long double Vr = entering->mat.FermiReal*1e-9; //Real component of the Fermi potential of entering material in [eV]
			long double kl = sqrt(2*neutron_mass_c2*Vr)/hbarc_squared; //limiting wave number defined in Heule's thesis after Eq 4.11
			long double transprob = 4*k1*k3/(k1 + k3)/(k1 + k3)*sqrt(1+(2*b*b*kl*kl)/(1+0.85*kl*w+2*kl*kl*w*w)); // transmission probability from Eq 4.26 of Heule's thesis
			reflprob = 1 - transprob;
		}
		//cout << " This is the reflection probability after MR modification = " << transprob << '\n';
		return reflprob;

	}

	void DoSnellTransmission(long double y1[6], long double y2[6], const long double normal[3], long double Estep) {
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2];
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
		long double k1 = sqrt(Enormal); // wavenumber in first solid 
		long double k2 = sqrt(Enormal - Estep); // wavenumber in second solid
		for (int i = 0; i < 3; i++)
			y2[i + 3] += (k2/k1 - 1)*(normal[i]*vnormal); // refract (scale normal velocity by k2/k1)
	}

	void DoAbsorption(long double x1, long double y1[6], long double &x2, long double y2[6], solid *entering) {
		ID = entering->ID; // set ID to absorbing solid
		x2 = x1; // set stopping point to track point right before collision
		for (int i = 0; i < 6; i++)
			y2[i] = y1[i];
	}

	void DoSpecularReflection(const long double normal[3], long double x1, long double y1[6], long double &x2, long double y2[6]) {
			
			//printf("Specular reflection! Erefl=%LG neV\n",Enormal*1e9);
			long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
			x2 = x1;
			for (int i = 0; i < 6; i++)
				y2[i] = y1[i];
			y2[3] -= 2*vnormal*normal[0]; // reflect velocity
			y2[4] -= 2*vnormal*normal[1];
			y2[5] -= 2*vnormal*normal[2];
	}
	
	void DoDiffuseReflection(long double x1, long double y1[6], long double &x2, long double y2[6], const long double normal[3]) {
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2];		
		long double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
		long double phi_r = mc->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
		long double theta_r = mc->SinCosDist(0, 0.5*pi);
		if (vnormal > 0) theta_r += pi; // if normal points out of volume rotate by 180 degrees
		x2 = x1;
		for (int i = 0; i < 6; i++)					y2[i] = y1[i];
		y2[3] = vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis
		y2[4] = vabs*sin(phi_r)*sin(theta_r);
		y2[5] = vabs*cos(theta_r);
		RotateVector(&y2[3],normal); // rotate coordinate system so that new z-axis lies on normal
//		printf("Diffuse reflection! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);

	}

	void DoDiffuseTransmission(long double x1, long double y1[6], long double &x2, long double y2[6], const long double normal[3], long double Estep) {
		//THIS FUNCTION NEEDS TO BE TESTED
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2];
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane
		long double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
		long double phi_r = mc->UniformDist(0, 2*pi); // generate random reflection angles (Lambert's law)
		long double theta_r = mc->SinCosDist(0, 0.5*pi);
		if (vnormal < 0) theta_r += pi; // if normal points into volume rotate by 180 degrees
		x2 = x1;
		for (int i = 0; i < 6; i++)					y2[i] = y1[i];
		long double k1 = sqrt(Enormal); // wavenumber in first solid 
		long double k2 = sqrt(Enormal - Estep); // wavenumber in second solid
		y2[3] = (k2/k1 - 1)*vabs*cos(phi_r)*sin(theta_r);	// new velocity with respect to z-axis, with scaling factor (k2/k1 -1)
		y2[4] = (k2/k1 - 1)*vabs*sin(phi_r)*sin(theta_r);
		y2[5] = (k2/k1 - 1)*vabs*cos(theta_r);
		RotateVector(&y2[3],normal); // rotate coordinate system so that new z-axis lies on normal
//		printf("Diffuse transmission! Erefl=%LG neV w_e=%LG w_s=%LG\n",Enormal*1e9,phi_r/conv,theta_r/conv);
	}

	void DoMicroRoughnessReflection(long double Estep, long double x1, long double y1[6], long double &x2, long double y2[6], const long double normal[3], solid *entering, MRHelperPointerless MRFunctions) {
		bool found = false;
		long double E = Ekin(&y1[3]); //Energy of Neutron
		long double theta_i = IncidentTheta(y1, normal);
		long double b = entering->mat.RMSRough*1e-9; // RMS roughness of entering material in [m]
		long double w = entering->mat.CorrelLength*1e-9; // Correlation length of entering material in [m]
		while (found == false) {
			long double phi_o_test = mc->UniformDist(-1*pi,pi); //this is the potential outgoing phi value for the transmitted neutron
			long double theta_o_test = mc->UniformDist(0,0.5*pi); //this is the potential outgoing theta value for the transmitted neutron
			long double prob_test = mc->UniformDist(0,0.5); //this is the probability value for the given phi and theta value above,
								     //which will be compared to see if it is above or below the actual value computed (0.5 used as max isn't available)
			double AngCut = entering->mat.AngCut;			
			long double prob_actual = MRFunctions.ProbIplus(E, Estep, theta_i, theta_o_test, phi_o_test, b, w, AngCut);
			if (prob_test <= prob_actual) {
				x2 = x1;
				for (int i = 0; i < 6; i++)
				y2[i] = y1[i];
				MicroRoughnessRotateVector(&y2[3], normal, theta_o_test-theta_i, pi+phi_o_test);	
				found = true;
				//cout << "This is the normal = (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << '\n';
				//cout << "This is the initial velocity vector = (" << y1[3] << ", " << y1[4] << ", " << y1[5] << ")" << '\n';
				//cout << "This is the MR reflected velocity vector = (" << y2[3] << ", " << y2[4] << ", " << y2[5] << ")" << '\n';
				//cout << "This neutron has been reflected according to the MR Model!" << '\n';
			}
		}

	}

	void DoMicroRoughnessTransmission(long double Estep, long double x1, long double y1[6], long double &x2, long double y2[6], const long double normal[3], solid *entering, MRHelperPointerless MRFunctions) {
		bool found = false;
		long double E = Ekin(&y1[3]); //Kinetic Energy of Neutron
		long double theta_i = IncidentTheta(y1,normal);
		long double b = entering->mat.RMSRough*1e-9; // RMS roughness of entering material in [m]
		long double w = entering->mat.CorrelLength*1e-9; // correlation length of entering material in [m]
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2];	
		long double Enormal = 0.5*m_n*vnormal*vnormal; // energy normal to reflection plane	
		while (found == false) {
			//long double vabs = sqrt(y1[3]*y1[3] + y1[4]*y1[4] + y1[5]*y1[5]);
			long double phi_o_test = mc->UniformDist(-1*pi,pi); //this is the potential outgoing phi value for the transmitted neutron
			long double theta_o_test = mc->UniformDist(0,0.5*pi); //this is the potential outgoing theta value for the transmitted neutron
			long double prob_test = mc->UniformDist(0,0.5); //this is the probability value for the given phi and theta value
									 //above which will be compared to see if it is above or below the actual value computed
			double AngCut = entering->mat.AngCut;
			long double prob_actual = MRFunctions.ProbIminus(E, Estep, theta_i, theta_o_test, phi_o_test, b, w, AngCut);
			if (prob_test <= prob_actual) {
				x2 = x1;
				for (int i = 0; i < 3; i++)
				y2[i] = y1[i];
				for (int i = 3; i < 6; i++) {//This is to scale the velocity according to the material it is travelling through
					long double k1 = sqrt(Enormal); // wavenumber in first solid (use only real part for transmission!)
					long double k2 = sqrt(Enormal - Estep); // wavenumber in second solid (use only real part for
					y2[i] = (k2/k1 - 1)*y2[i];
				}	
				MicroRoughnessRotateVector(&y2[3], normal,pi-theta_o_test-theta_i, pi+phi_o_test);
				found = true;
				cout << "This neutron has been transmitted according to the MR Model!" << '\n';
			}
		}	

	}

	long double IncidentTheta(long double y1[6],const long double normal[3]) {
		long double vnormal = y1[3]*normal[0] + y1[4]*normal[1] + y1[5]*normal[2]; // velocity normal to reflection plane
		long double vmagnitude = sqrt(y1[3]*y1[3]+y1[4]*y1[4]+y1[5]*y1[5]); //the magnitude of the velocity
		long double normalmagnitude = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]); //the magnitude of the normal vector
		long double theta_i = pi - acos(vnormal/(vmagnitude*normalmagnitude)); //the incident neutron angle in radians
		return theta_i;
	}


};


#endif // NEUTRON_H_
