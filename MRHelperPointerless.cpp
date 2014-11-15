//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
////
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// This file contains the source code of various functions all related to the
// calculation of microroughness.
//
// see A. Steyerl, Z. Physik 254 (1972) 169.
//
// A description of the functions can be found in the corresponding header file
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MRHelperPointerless.h"
#include <complex>
#include <math.h>
#include <stdlib.h>
#include "globals.h"
#include <stdio.h>
#include <iostream>

using namespace std;

MRHelperPointerless* MRHelperPointerless::fpInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor

MRHelperPointerless::MRHelperPointerless() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MRHelperPointerless::~MRHelperPointerless()
{double IntIplus(double, double, double, int, int, double, double, double);
  delete fpInstance;
  fpInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MRHelperPointerless* MRHelperPointerless::GetInstance()
{
  if (fpInstance == 0) fpInstance = new MRHelperPointerless;
  return fpInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::S2(double costheta2, double klk2)
{
  // case 1: radicand is positive,
  // case 2: radicand is negative, cf. p. 174 of the Steyerl paper

  if (costheta2>=klk2)
     return 4*costheta2/(2*costheta2-klk2+2*sqrt(costheta2*(costheta2-klk2)));
  else
     return 4*costheta2/klk2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::SS2(double costhetas2, double klks2)
{
  // cf. p. 175 of the Steyerl paper

  return 4*costhetas2/
                   (2*costhetas2+klks2+2*sqrt(costhetas2*(costhetas2+klks2)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::Fmu(double k2, double thetai,
                                        double thetao, double phio,
                                        double b2, double w2,
                                        double AngCut)
{
  double mu_squared;

  // Checks if the distribution is peaked around the specular direction

  if ((abs(thetai-thetao)<AngCut) && (abs(phio)<AngCut))
    mu_squared=0.;
  else
    {
      // cf. p. 177 of the Steyerl paper

      double sinthetai=sin(thetai);
      double sinthetao=sin(thetao);
      mu_squared=k2*
          (sinthetai*sinthetai+sinthetao*sinthetao-
                   2.*sinthetai*sinthetao*cos(phio));
    }

  // cf. p. 177 of the Steyerl paper

  return b2*w2/twopi*exp(-mu_squared*w2/2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::FmuS(double k, double kS,
                                         double thetai, double thetaSo,
                                         double phiSo,
                                         double b2, double w2,
                                         double AngCut, double thetarefract)
{
  double mu_squared;

  // Checks if the distribution is peaked around the direction of
  // unperturbed refraction

  if ((abs(thetarefract-thetaSo)<AngCut) && (abs(phiSo)<AngCut))
    mu_squared=0.;
  else
    {
      double sinthetai=sin(thetai);
      double sinthetaSo=sin(thetaSo);

      // cf. p. 177 of the Steyerl paper
      mu_squared=k*k*sinthetai*sinthetai+kS*kS*sinthetaSo*sinthetaSo-
                 2.*k*kS*sinthetai*sinthetaSo*cos(phiSo);
    }

  // cf. p. 177 of the Steyerl paper

  return b2*w2/twopi*exp(-mu_squared*w2/2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::IntIplus(double E, double fermipot,
                                             double theta_i, int AngNoTheta,
                                             int AngNoPhi, double b2,
                                             double w2,
                                             double AngCut)
{
  

  // max_theta_o saves the theta-position of the max probability,
  // the previous value is saved in a_max_theta_o

  double a_max_theta_o, max_theta_o=theta_i, a_max_phi_o, max_phi_o=0.;

  // max_phi_o saves the phi-position of the max probability,
  // the previous value is saved in a_max_phi_o

  // Definition of the stepsizes in theta_o and phi_o

  double theta_o;
  double phi_o;
  double Intens;
  double ang_steptheta=90.*conv/(AngNoTheta-1);
  double ang_stepphi=360.*conv/(AngNoPhi-1);
  double costheta_i=cos(theta_i);
  double costheta_i_squared=costheta_i*costheta_i;

  // (k_l/k)^2
  double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                                 hbarc_squared*fermipot*fermipot;

  // (k_l/k)^2
  double klk2=fermipot/E;

  double costheta_o_squared;

  // k^2
  double k2=2*neutron_mass_c2*E/hbarc_squared;

  double wkeit=0.;

  // Loop through theta_o
 
  for (theta_o=0.*conv; theta_o<=90.*conv+1e-6; theta_o+=ang_steptheta)
    {
      costheta_o_squared=cos(theta_o)*cos(theta_o);

      // Loop through phi_o

      for (phi_o=-180.*conv; phi_o<=180.*conv+1e-6; phi_o+=ang_stepphi)
	{
          //calculates probability for a certain theta_o,phi_o pair

	  Intens=kl4d4/costheta_i*S2(costheta_i_squared,klk2)*
                 S2(costheta_o_squared,klk2)*
                 Fmu(k2,theta_i,theta_o,phi_o,b2,w2,AngCut)*sin(theta_o);

          //cout << "S2:  " << S2(costheta_i_squared,klk2) << endl;
          //cout << "S2:  " << S2(costheta_o_squared,klk2) << endl;
          //cout << "Fmu: " << Fmu(k2,theta_i,theta_o,phi_o,b2,w2,AngCut)*sin(theta_o) << endl;
          // checks if the new probability is larger than the
          // maximum found so far

          
          // Adds the newly found probability to the integral probability

	  wkeit+=Intens*ang_steptheta*ang_stepphi;
	}
    }

  // Fine-Iteration to find maximum of distribution
  // only if the energy is not zero

  if (E>1e-10*eV)
    {

  // Break-condition for refining

  while ((ang_stepphi>=AngCut*AngCut) || (ang_steptheta>=AngCut*AngCut))
    {
      a_max_theta_o=max_theta_o;
      a_max_phi_o=max_phi_o;
      ang_stepphi /= 2.;
      ang_steptheta /= 2.;

      //cout << ang_stepphi/conv << ", "
      //       << ang_steptheta/conv << ","
      //       << AngCut/conv << endl;

      for (theta_o=a_max_theta_o-ang_steptheta;
           theta_o<=a_max_theta_o-ang_steptheta+1e-6;
           theta_o+=ang_steptheta)
	{
	  //cout << "theta_o: " << theta_o/conv << endl;
	  costheta_o_squared=cos(theta_o)*cos(theta_o);
	  for (phi_o=a_max_phi_o-ang_stepphi;
               phi_o<=a_max_phi_o+ang_stepphi+1e-6;
               phi_o+=ang_stepphi)
	    {
	      //cout << "phi_o: " << phi_o/conv << endl;
	      Intens=kl4d4/costheta_i*S2(costheta_i_squared, klk2)*
                     S2(costheta_o_squared,klk2)*
                     Fmu(k2,theta_i,theta_o,phi_o,b2,w2,AngCut)*sin(theta_o);
	      
	    }
	}
    }
    }
  return wkeit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::IntIminus(double E, double fermipot,
                                              double theta_i,
                                              int AngNoTheta,
                                              int AngNoPhi, double b2,
                                              double w2, 
                                              double AngCut)
{
  double a_max_thetas_o, max_thetas_o = theta_i;
  double a_max_phis_o, max_phis_o = 0.;
  double thetas_o;
  double phis_o;
  double IntensS;
  double ang_steptheta=180.*conv/(AngNoTheta-1);
  double ang_stepphi=180.*conv/(AngNoPhi-1);
  double costheta_i=cos(theta_i);
  double costheta_i_squared=costheta_i*costheta_i;


  double wkeit=0.;

  if (E*costheta_i_squared < fermipot) return wkeit; //added cos^2(theta_i) in front of E, to make it energy normal to surface

  //k_l^4/4
  double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                 hbarc_squared*fermipot*fermipot;
  // (k_l/k)^2
  double klk2=fermipot/E;

  // (k_l/k')^2
  double klks2=fermipot/(E-fermipot);

  // k'/k
  double ksdk=sqrt((E-fermipot)/E);

  double costhetas_o_squared;

  // k
  double k=sqrt(2*neutron_mass_c2*E/hbarc_squared);

  // k'
  double kS=ksdk*k;

  for (thetas_o=0.*conv; thetas_o<=90.*conv+1e-6; thetas_o+=ang_steptheta)
    {
      costhetas_o_squared=cos(thetas_o)*cos(thetas_o);

      for (phis_o=-180.*conv; phis_o<=180.*conv+1e-6; phis_o+=ang_stepphi)
	{
          //cf. Steyerl-paper p. 176, eq. 21
	  if (costhetas_o_squared>=-klks2) {

            //calculates probability for a certain theta'_o, phi'_o pair

            double thetarefract = thetas_o;
            if (std::abs(sin(theta_i)/ksdk) <= 1.)
                                        thetarefract = asin(sin(theta_i)/ksdk);

	    IntensS = kl4d4/costheta_i*ksdk*S2(costheta_i_squared, klk2)*
                      SS2(costhetas_o_squared,klks2)*
                      FmuS(k,kS,theta_i,thetas_o,phis_o,b2,w2,AngCut,thetarefract)*
                      sin(thetas_o);
	  } else {
	    IntensS=0.;
          }
            // checks if the new probability is larger than
            // the maximum found so far
	  
	  wkeit+=IntensS*ang_steptheta*ang_stepphi;
	}
    }

  // Fine-Iteration to find maximum of distribution

  if (E>1e-10*eV)
    {

  // Break-condition for refining

  while (ang_stepphi>=AngCut*AngCut || ang_steptheta>=AngCut*AngCut)
    {
      a_max_thetas_o=max_thetas_o;
      a_max_phis_o=max_phis_o;
      ang_stepphi /= 2.;
      ang_steptheta /= 2.;
      //cout << ang_stepphi/conv << ", " << ang_steptheta/conv 
      //       << ", " << AngCut/conv << endl;
      for (thetas_o=a_max_thetas_o-ang_steptheta;
           thetas_o<=a_max_thetas_o-ang_steptheta+1e-6;
           thetas_o+=ang_steptheta)
	{
	  costhetas_o_squared=cos(thetas_o)*cos(thetas_o);
	  for (phis_o=a_max_phis_o-ang_stepphi;
               phis_o<=a_max_phis_o+ang_stepphi+1e-6;
               phis_o+=ang_stepphi)
	    {
              double thetarefract = thetas_o;
              if (std::abs(sin(theta_i)/ksdk) <= 1.)
                                       thetarefract = asin(sin(theta_i)/ksdk);

	      IntensS=kl4d4/costheta_i*ksdk*S2(costheta_i_squared, klk2)*
                      SS2(costhetas_o_squared,klks2)*
                      FmuS(k,kS,theta_i,thetas_o,phis_o,b2,w2,AngCut,thetarefract)*
                      sin(thetas_o);
	   
	    }
	}
    }
    }
  return wkeit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::ProbIplus(double E, double fermipot,
                                              double theta_i,
                                              double theta_o,
                                              double phi_o,
                                              double b, double w,
                                              double AngCut)
{
  //k_l^4/4
  double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                 hbarc_squared*fermipot*fermipot;

  // (k_l/k)^2
  double klk2=fermipot/E;

  double costheta_i=cos(theta_i); 
  double costheta_o=cos(theta_o);

  // eq. 20 on page 176 in the steyerl paper

  return kl4d4/costheta_i*S2(costheta_i*costheta_i, klk2)*
         S2(costheta_o*costheta_o,klk2)*
   Fmu(2*neutron_mass_c2*E/hbarc_squared,theta_i,theta_o,phi_o,b*b,w*w,AngCut)*
   sin(theta_o);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double MRHelperPointerless::ProbIminus(double E, double fermipot,
                                               double theta_i,
                                               double thetas_o,
                                               double phis_o, double b,
                                               double w, double AngCut)
{
  //k_l^4/4
  double kl4d4=neutron_mass_c2/hbarc_squared*neutron_mass_c2/
                 hbarc_squared*fermipot*fermipot;
  // (k_l/k)^2
  double klk2=fermipot/E;

  // (k_l/k')^2
  double klks2=fermipot/(E-fermipot);

  if (E < fermipot) {
     cout << " ProbIminus E < fermipot " << endl;
     return 0.;
  }

  // k'/k
  double ksdk=sqrt((E-fermipot)/E);

  // k
  double k=sqrt(2*neutron_mass_c2*E/hbarc_squared);

  // k'/k
  double kS=ksdk*k;

  double costheta_i=cos(theta_i);
  double costhetas_o=cos(thetas_o);

  // eq. 20 on page 176 in the steyerl paper

  double thetarefract = thetas_o;
  if(std::abs(sin(theta_i)/ksdk) <= 1.)thetarefract = asin(sin(theta_i)/ksdk);

  return kl4d4/costheta_i*ksdk*S2(costheta_i*costheta_i, klk2)*
         SS2(costhetas_o*costhetas_o,klks2)*
         FmuS(k,kS,theta_i,thetas_o,phis_o,b*b,w*w,AngCut,thetarefract)*
         sin(thetas_o);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
