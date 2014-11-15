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
//---------------------------------------------------------------------------
//
// ClassName:   MicroRoughnessHelper
//
// Class description:
//
// This file contains the headers of various functions all related to the
// calculation of microroughness.
// see A. Steyerl, Z. Physik 254 (1972) 169.
//
// Angular cut: for angles which are closer to the specular direction than a
// certain value (0.01°), the probability is set to 0 in order to avoid a
// hang-up at the generation of the polar angle due to a very sharp angular
// distribution
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 12-05-14, adopted from Stefan Heule (PSI) Thesis by P.Gumplinger
//           http://.web.psi.ch/papers/stefanheule_thesis2008.pdf
//           reported in F. Atchison et al., Eur. Phys. J. A 44, 23–29 (2010)
//           Thanks to Geza Zsigmond

#ifndef MICROROUGHNESSHELPER_H
#define MICROROUGHNESSHELPER_H 1



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MicroRoughnessHelper
{
  public:  // with description

    static MicroRoughnessHelper* GetInstance();

    MicroRoughnessHelper();
    ~MicroRoughnessHelper();



// Transmitted intensity with k-vector in vacuum

// arguments:
//         1) cos(theta)^2,
//         2) (k_l/k)^2

    double S2(double, double);

// Transmitted intensity with k-vector within the medium

// arguments:
//         1) cos(theta')^2,
//         2) (k_l/k')^2

    double SS2(double, double);

// Fourier-tranform of the autocorrelation function with k-vector in vacuum

// arguments:
//         1) k^2,
//         2) theta_i,
//         3) theta_o,
//         4) phi_o,
//         5) b^2,
//         6) w^2,
//         7) angular cut

    double Fmu(double, double, double,
                 double, double, double, double);

// Fourier-tranform of the autocorrelation function with k-vector within
// the medium

// arguments:
//         1) k,
//         2) k',
//         3) theta_i,
//         4) theta'_o,
//         5) phi'_o,
//         6) b^2,
//         7) w^2,
//         8) angular cut
//         9) theta_refract

    double FmuS(double, double, double, double,
                  double, double, double, double, double);

// Integral probability for non-specular reflection

// arguments:
//         1) E,
//         2) V_F,
//         3) theta_i,
//         4) number of angles theta_o for which the probability is calculated,
//         5) number of angles phi_o for which the probability is calculated,
//         6) b^2,
//         7) w^2,
//         8) pointer to double array with max values of the probability,
//         9) angular cut

    double IntIplus(double, double, double, int, int,
                      double, double, double*, double);

// Probability of non-specular reflection with the microroughness model

// arguments:
//         1) E,
//         2) V_F,
//         3) theta_i,
//         4) theta_o,
//         5) phi_o,
//         6) b,
//         7) w,
//         8) angular cut

    double ProbIplus (double, double, double, double,
                        double, double, double, double);

// Integral probability for non-specular transmission

// arguments:
//         1) E,
//         2) V_F,
//         3) theta_i,
//         4) number of angles theta_o for which the probability is calculated,
//         5) number of angles phi_o for which the probability is calculated,
//         6) b^2,
//         7) w^2,
//         8) pointer to double array with max values of the probability,
//         9) angular cut

    double IntIminus(double, double, double, int, int,
                       double, double, double*, double);

// Probability of non-specular transmission with the microroughness model

// arguments:
//         1) E,
//         2) V_F,
//         3) theta_i,
//         4) theta'_o,
//         5) phi'_o,
//         6) b,
//         7) w,
//         8) angular cut

    double ProbIminus (double, double, double, double,
                         double, double, double, double);

  private:

    static MicroRoughnessHelper* fpInstance;

};

#endif // MICROROUGHNESSHELPER_H
