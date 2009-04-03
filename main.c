#include "main.h"
//#include "lines.h"
/*
b = 0  kein bereich alle zeilen werden ignoriert bis auf [ startende
b = 1  globale settings
b = 2  ramp up
b = 3  full field
b = 4  ramp down
*/

void ConfigInit(void){
	FILE *cfg = NULL;
	int n=0,b=0;
	char line[1024];
	
	string path;
	path = inpath + "/config.in";
	cfg = fopen(path.c_str(),mode_r);
	if(cfg == NULL){  
		exit(-1);
		return;
	}
	
	/* setting default values */
	protneut = NEUTRON;
	runge = 2;
	polarisation = 1;
	spinflipcheck = 0;
	BruteForce = 0;
	neutdist = 0;
	wanddicke = 0.01;
	reflekt = 1;
	diffuse = 3;
	bfeldwahl = 0;
	ausgabewunsch = 2;
	ausgabewunschsave = 2;
	MonteCarlo = 1;
	MonteCarloAnzahl = 1;
	nodelay = 1;
	M = m_n;
	mu_n = mu_nSI / ele_e * -1;
	mumB = mu_n/M;
	decay = 0;
	/*end default values*/
	
	// we want do find some keywords in the config file bah[][] contains the possible variables and bah2[][] the regions
	char bah[28][17]={{'s','l','i','t'},{'B','r','u','t','e','F','o','r','c','e'},{'r','e','f','l','e','k','t'},
		{'s','p','i','n','f','l','i','p','c','h','e','c','k'},{'p','r','o','t','n','e','u','t'},{'r','u','n','g','e'},{'p','o','l','a','r','i','s','a','t','i','o','n'},
		{'n','e','u','t','d','i','s','t'},{'w','a','n','d','d','i','c','k','e'},{'d','i','f','f','u','s','e'},{'b','f','e','l','d','w','a','h','l'},
		{'a','u','s','g','a','b','e','w','u','n','s','c','h'},{'M','o','n','t','e','C','a','r','l','o',' '},{'M','o','n','t','e','C','a','r','l','o','A','n','z','a','h','l'},
		{'R','a','m','p','U','p','T','i','m','e'},{'R','a','m','p','D','o','w','n','T','i','m','e'},{'F','u','l','l','F','i','e','l','d','T','i','m','e'},{'C','l','e','a','n','i','n','g','T','i','m','e'},{'E','m','p','t','y','i','n','g','T','i','m','e'},{'s','t','o','r','a','g','e','t','i','m','e'},
		{'B','F','T','a','r','g','e','t','B'},{'d','e','c','a','y'},{'a','b','s','o','r','b','e','r'},{'F','i','l','l','i','n','g','T','i','m','e'},{'D','e','t','O','p','e','n'},{'B','F','e','l','d','S','k','a','l','G','l','o','b','a','l'},{'E','F','e','l','d','S','k','a','l'},
			{'e','x','p','m','o','d','e'}};
			
	char bah2[7][10]={{'g','l','o','b','a','l'},{'r','a','m','p','u','p'},{'f','u','l','l','f','i','e','l','d'},{'r','a','m','p','d','o','w','n'},{'f','i','l','l','i','n','g'},{'c','o','u','n','t','i','n','g'},{'c','l','e','a','n','i','n','g'}};
	
	do // read the file
	{
		n=0;
		fgets(line,1024,cfg); // get a line
		if(line[0]=='#') continue; // if its a comment goto next line
		if((b == 0) && (line[0]!='[')) continue; // if no region specified and this line doesn't start with [ goto next line
		if(line[0]=='['){ // no region specified but we found [ as a first char now we go to find out which region 
			if(strstr(line,bah2[0]) != NULL){
				b=1;
				continue;
			}else if(strstr(line,bah2[1]) != NULL){
				b=2;
				continue;
			}else if(strstr(line,bah2[2]) != NULL){
				b=3;
				continue;
			}else if(strstr(line,bah2[3]) != NULL){
				b=4;
				continue;
			}else if(strstr(line,bah2[4]) != NULL){
				b=5;
				continue;
			}else if(strstr(line,bah2[5]) != NULL){
				b=6;
				continue;
				}else if(strstr(line,bah2[6]) != NULL){
				b=7;
				continue;
			}else{
				b=0;
				continue;
			}
		}else{ // we have a region and want to analyse the line
			switch(b){
				case 1:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %hu",&slit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %hu",&BruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&reflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&spinflipcheck);
					else if(strstr(line,bah[4]) != NULL)
						sscanf(line,"%*s %i",&protneut);
					else if(strstr(line,bah[5]) != NULL)
						sscanf(line,"%*s %d",&runge);
					else if(strstr(line,bah[6]) != NULL)
						sscanf(line,"%*s %i",&polarisation);
					else if(strstr(line,bah[7]) != NULL)
						sscanf(line,"%*s %i",&neutdist);
					else if(strstr(line,bah[8]) != NULL)
						sscanf(line,"%*s %5LG",&wanddicke);
					else if(strstr(line,bah[9]) != NULL)
						sscanf(line,"%*s %i",&diffuse);
					else if(strstr(line,bah[10]) != NULL)
						sscanf(line,"%*s %i",&bfeldwahl);
					else if(strstr(line,bah[11]) != NULL)
						sscanf(line,"%*s %i",&ausgabewunsch);
					else if(strstr(line,bah[12]) != NULL)
						sscanf(line,"%*s %i",&MonteCarlo);
					else if(strstr(line,bah[13]) != NULL)
						sscanf(line,"%*s %i",&MonteCarloAnzahl);
					else if(strstr(line,bah[14]) != NULL)
						sscanf(line,"%*s %7LG",&RampUpTime);
					else if(strstr(line,bah[15]) != NULL)
						sscanf(line,"%*s %7LG",&RampDownTime);
					else if(strstr(line,bah[16]) != NULL)
						sscanf(line,"%*s %7LG",&FullFieldTime);
					else if(strstr(line,bah[17]) != NULL)
						sscanf(line,"%*s %7LG",&CleaningTime);
					else if(strstr(line,bah[18]) != NULL)
						sscanf(line,"%*s %7LG",&EmptyingTime);
		            else if(strstr(line,bah[19]) != NULL)
						sscanf(line,"%*s %7LG",&storagetime);
					else if(strstr(line,bah[20]) != NULL)
						sscanf(line,"%*s %3LG",&BFTargetB);
					else if(strstr(line,bah[21]) != NULL)
						sscanf(line,"%*s %hu",&decay);
					else if(strstr(line,bah[22]) != NULL)
						sscanf(line,"%*s %d",&AbsorberChoice);
					else if(strstr(line,bah[23]) != NULL)
						sscanf(line,"%*s %7LG",&FillingTime);
					else if(strstr(line,bah[25]) != NULL)
						sscanf(line,"%*s %7LG",&BFeldSkalGlobal);
					else if(strstr(line,bah[26]) != NULL)
						sscanf(line,"%*s %7LG",&EFeldSkal);
					else if(strstr(line,bah[27]) != NULL)
						sscanf(line,"%*s %d",&expmode);
				break;
				
				case 2:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %i",&ruslit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %i",&ruBruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&rureflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&ruspinflipcheck);
					else if(strstr(line,bah[24]) != NULL)
						sscanf(line,"%*s %i",&ruDetOpen);
				break;
			
				case 3:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %i",&ffslit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %i",&ffBruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&ffreflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&ffspinflipcheck);
					else if(strstr(line,bah[24]) != NULL)
						sscanf(line,"%*s %i",&ffDetOpen);
				break;
				
				case 4:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %i",&rdslit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %i",&rdBruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&rdreflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&rdspinflipcheck);
					else if(strstr(line,bah[24]) != NULL)
						sscanf(line,"%*s %i",&rdDetOpen);
				break;
					
				case 5:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %i",&fislit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %i",&fiBruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&fireflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&fispinflipcheck);
					else if(strstr(line,bah[24]) != NULL)
						sscanf(line,"%*s %i",&fiDetOpen);
				break;
					
				case 6:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %i",&coslit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %i",&coBruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&coreflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&cospinflipcheck);
					else if(strstr(line,bah[24]) != NULL)
						sscanf(line,"%*s %i",&coDetOpen);
				break;
					
					case 7:
					if(strstr(line,bah[0]) != NULL)
						sscanf(line,"%*s %i",&clslit);
					else if(strstr(line,bah[1]) != NULL)
						sscanf(line,"%*s %i",&clBruteForce);
					else if(strstr(line,bah[2]) != NULL)
						sscanf(line,"%*s %i",&clreflekt);
					else if(strstr(line,bah[3]) != NULL)
						sscanf(line,"%*s %i",&clspinflipcheck);
					else if(strstr(line,bah[24]) != NULL)
						sscanf(line,"%*s %i",&clDetOpen);
				break;
				
				case 0:
					continue;
				break;
				
				default:
					continue;
				break;
			}
		}
		
	}while(!feof(cfg));
	
	if(protneut == NEUTRON)
	{
		M=m_n;   // [eV/c^2]
		// make absorber choice
		if(AbsorberChoice==1) //PE
		{ Mf=FPrealPE; Pf=FPimPE;	}
		else if(AbsorberChoice==2) //Ti
		{ Mf=FPrealTi; Pf=FPimTi;	} 
		if(polarisation==POLARISATION_GOOD){
			hfs = -1;
			mu_n=hfs * mu_nSI / ele_e; //=> mu_n is negative
			mumB= mu_n/M;  // [c^2/T]
		}else if(polarisation==POLARISATION_BAD){  // in ev/T
			hfs = 1;
			mu_n=hfs * mu_nSI / ele_e; //=> mu_n is positive
			mumB= mu_n/M;  // [c^2/T]
		} else if (polarisation==POLARISATION_NONE){
			hfs = 0;
			mu_n=0;
			mumB= mu_n/M;  // [c^2/T] 
		}  // in ev/T
		if (neutdist == 1) prepndist(1);
	}
	else if(protneut==PROTON)
	{
		M=m_p;
		Qm0=1.0/M;
	}
	else if (protneut == ELECTRONS)
	{
		M=m_e;
		Qm0=-1.0/M;
	}
	
	if((BruteForce==1)||(ausgabewunsch==1)||(neutdist==1))
	{
		SaveIntermediate=1;
		kmax=KMDEF;
	}
	
	Efeldwahl=0;
	ausgabewunschsave=ausgabewunsch;
	BFeldSkalGlobalSave=BFeldSkalGlobal;
	EFeldSkalSave=EFeldSkal;
	fclose(cfg);
	return;
}



//======== initalizes initial values from record to used variables ==============================================
void initialStartbed()
{	switch(protneut)
	{	case NEUTRON:	EnergieS = nini.EnergieS;
						dEnergie = nini.dEnergie;
						EnergieE = nini.EnergieE;
						z_ns = nini.zs;
						dz_n = nini.dz;
						z_ne = nini.ze;
						r_ns = nini.rs;
						dr_n = nini.dr;
						r_ne = nini.re;						
						phis = nini.phis;
						dphi = nini.dphi;
						phie = nini.phie;
						alphas = nini.alphas;
						dalpha = nini.dalpha;
						alphae = nini.alphae;
						gammas = nini.gammas;
						dgamma = nini.dgamma;
						gammae = nini.gammae;
						delx = nini.delx;
						xend = nini.xend;
						break;	
		case PROTON:	EnergieS = pini.EnergieS;
						dEnergie = pini.dEnergie;
						EnergieE = pini.EnergieE;
						z_ns = pini.zs;
						dz_n = pini.dz;
						z_ne = pini.ze;
						r_ns = pini.rs;
						dr_n = pini.dr;
						r_ne = pini.re;						
						phis = pini.phis;
						dphi = pini.dphi;
						phie = pini.phie;
						alphas = pini.alphas;
						dalpha = pini.dalpha;
						alphae = pini.alphae;
						gammas = pini.gammas;
						dgamma = pini.dgamma;
						gammae = pini.gammae;
						delx = pini.delx;
						xend = pini.xend;
						break;
		case ELECTRONS:	EnergieS = eini.EnergieS;
						dEnergie = eini.dEnergie;
						EnergieE = eini.EnergieE;
						z_ns = eini.zs;
						dz_n = eini.dz;
						z_ne = eini.ze;
						r_ns = eini.rs;
						dr_n = eini.dr;
						r_ne = eini.re;
						phis = eini.phis;
						dphi = eini.dphi;
						phie = eini.phie;
						alphas = eini.alphas;
						dalpha = eini.dalpha;
						alphae = eini.alphae;
						gammas = eini.gammas;
						dgamma = eini.dgamma;
						gammae = eini.gammae;
						delx = eini.delx;
						xend = eini.xend;
						break;
	}
	//return;
}
//======== end of innitialStartbed ==============================================================================

//======== read in of the initial values (out of all3inone.in) and the dimensions (out of dimensions.in) ========
void Startbed(int k)
{	int i = 0, ncont = 42;
    char cline[200];    
	string path;

	// setting path for all3inone.in
	if(protneut==PROTON || protneut==NEUTRON || protneut==ELECTRONS)
	{	path = inpath + "/all3inone.in";
	}
    else
    {	exit(-1);
    }

	// creating 'inistream' to all3inone.in
	FILE *inistream = fopen(path.c_str(), mode_r);
	if(inistream == NULL)
	{	cout << "Can't open " << path << "\n";
		exit(-1);
	}
	
	// looping  over the lines of all3inone.in and reading them in
	for(i = 1; i < 27; i++)
	{	fgets(cline, 200, inistream);
		switch (i)
		{	case  2:	ncont = sscanf(cline, "%LG %LG %LG ", &DiceRodField, &RodFieldMultiplicator, &Ibar);
						break;
			case  4:	ncont = sscanf(cline, "%LG %LG %LG ", &nini.EnergieS, &nini.dEnergie, &nini.EnergieE);
						break;
			case  5:	ncont = sscanf(cline, "%LG %LG %LG ", &nini.zs, &nini.dz, &nini.ze);
						break;
			case  6:	ncont = sscanf(cline, "%LG %LG %LG ", &nini.rs, &nini.dr , &nini.re);
						break;
			case  7:	ncont = sscanf(cline, "%LG %LG %LG ", &nini.phis, &nini.dphi, &nini.phie);
						break;			
			case  8:	ncont = sscanf(cline, "%LG %LG %LG ", &nini.alphas, &nini.dalpha, &nini.alphae);
						break;
			case  9:	ncont = sscanf(cline, "%LG %LG %LG ", &nini.gammas, &nini.dgamma, &nini.gammae);
						break;
			case 10:	ncont = sscanf(cline, "%LG %LG ", &nini.delx, &nini.xend);
						break;
			case 12:	ncont = sscanf(cline, "%LG %LG %LG ", &pini.EnergieS, &pini.dEnergie, &pini.EnergieE);
						break;
			case 13:	ncont = sscanf(cline, "%LG %LG %LG ", &pini.zs, &pini.dz, &pini.ze);
						break;
			case 14:	ncont = sscanf(cline, "%LG %LG %LG ", &pini.rs, &pini.dr , &pini.re);
						break;
			case 15:	ncont = sscanf(cline, "%LG %LG %LG ", &pini.phis, &pini.dphi, &pini.phie);
						break;			
			case 16:	ncont = sscanf(cline, "%LG %LG %LG ", &pini.alphas, &pini.dalpha, &pini.alphae);
						break;
			case 17:	ncont = sscanf(cline, "%LG %LG %LG ", &pini.gammas, &pini.dgamma, &pini.gammae);
						break;
			case 18:	ncont = sscanf(cline, "%LG %LG ", &pini.delx, &pini.xend);
						break;
			case 20:	ncont = sscanf(cline, "%LG %LG %LG ", &eini.EnergieS, &eini.dEnergie, &eini.EnergieE);
						break;
			case 21:	ncont = sscanf(cline, "%LG %LG %LG ", &eini.zs, &eini.dz, &eini.ze);
						break;
			case 22:	ncont = sscanf(cline, "%LG %LG %LG ", &eini.rs, &eini.dr , &eini.re);
						break;
			case 23:	ncont = sscanf(cline, "%LG %LG %LG ", &eini.phis, &eini.dphi, &eini.phie);
						break;			
			case 24:	ncont = sscanf(cline, "%LG %LG %LG ", &eini.alphas, &eini.dalpha, &eini.alphae);
						break;
			case 25:	ncont = sscanf(cline, "%LG %LG %LG ", &eini.gammas, &eini.dgamma, &eini.gammae);
						break;
			case 26:	ncont = sscanf(cline, "%LG %LG ", &eini.delx, &eini.xend);
						break;
		}
		if(ncont < 1) printf("an error occourd while reading the %i. item in line %i of all3inone.in", ncont, i);
    }
    
    // closing 'inistream' to all3inone.in
 	fclose(inistream);
	
	// setting path for dimensions.in
	path = inpath + "/dimensions.in";
	
	// creating 'dimstream' to dimensions.in	
	FILE *dimstream = fopen (path.c_str(), mode_r);
	if (dimstream == NULL)
	{	cout << "Can't open " << path << "\n";
		exit(-1);
	}   
	
	// looping  over the lines of dimensions.in and reading them in
	for(i = 1; i < 29; i++)
	{	fgets(cline, 200, dimstream);
		switch(i)
		{	case  2:	ncont = sscanf(cline, "%LG %LG ", &zmin, &zmax);
						break;
			case  3:	ncont = sscanf(cline, "%LG %LG ", &rmin, &rmax);
						break;
			case  4:	ncont = sscanf(cline, "%LG ", &detz);
						break;
			case  5:	ncont = sscanf(cline, "%LG %LG ", &detrmin, &detrmax);
						break;
			case  6:	ncont = sscanf(cline, "%LG %LG ", &abszmin, &abszmax);
						break;
			case  7:	ncont = sscanf(cline, "%LG %LG ", &absrmin, &absrmax);
						break;
			case  8:	ncont = sscanf(cline, "%LG %LG ", &absphimin, &absphimax);
						break;
			case  9:	ncont = sscanf(cline, "%LG ", &hlid);	
						break;
			case 11:	ncont = sscanf(cline, "%LG %LG", &FPrealNocado, &FPimNocado);
						break;
			case 12:	ncont = sscanf(cline, "%LG %LG", &FPrealPE, &FPimPE);
						break;
			case 13:	ncont = sscanf(cline, "%LG %LG", &FPrealTi, &FPimTi);
						break;
			case 14:	ncont = sscanf(cline, "%LG %LG", &FPrealCu, &FPimCu);
						break;
			case 15:	ncont = sscanf(cline, "%LG %LG", &FPrealCsI, &FPimCsI);
						break;
			case 16:	ncont = sscanf(cline, "%LG %LG", &FPrealDLC, &FPimDLC);
						break;						
		}
		if(ncont < 1) printf("an error occourd while reading the %i. item in line %i of dimensions.in", ncont, i);				
	}
	
	// closing 'dimstream' to all3inone.in
	fclose(dimstream);
	
	// initalizes initial values from record to used variables
	initialStartbed();
	
	// the inner wall thickness is set to the same size then the outer, but if there is no inner wall (/coils) the thickness is set to zero.
	// (mainly used for AbEx)
	wandinnen = wanddicke;
	if(rmin <= 0.0) wandinnen = 0.0;

	// logging dimensions to *log.out
	fprintf(LOGSCR, "dimensions.in: rmin = %.17LG, rmax = %.17LG, zmin = %.17LG, zmax = %.17LG\n"
	                "detz = %.17LG, detrmin = %.17LG, detrmax = %.17LG, hlid = %.17LG",
	                rmin, rmax, zmin, zmax,
	                detz, detrmin, detrmax, hlid);

	// writing parameters to screen
	printf("\nStart parameters:\n"
	       "Energy (min, step, max): %17LG neV/eV, %.17LG neV/eV, %.17LG neV/eV\n"
	       "Maximum runtime: %.17LG s\n"
	       "r (min, step, max): %.17LG m, %.17LG m, %.17LG m\n"
	       "z (min, step, max): %.17LG m, %.17LG m, %.17LG m\n"
	       "phi (min, step, max): %.17LG m, %.17LG m, %.17LG m\n"
	       "Alpha (min, step, max): %.17LG Grad, %.17LG Grad, %.17LG Grad\n"
	       "Gamma (min, step, max): %.17LG Grad, %.17LG Grad, %.17LG Grad\n"
	       "Filling time: %LG s\n"
	       "Cleaning time: %.17LG s\n"
	       "Ramp up time: %.17LG s\n"
	       "Full field time: %.17LG s\n"
	       "RampDownTime: %.17LG s\n"
	       "B field scaling factor: %.17LG s\n"
	       "E field scaling factor: %.17LG s\n",
	       EnergieS, dEnergie, EnergieE,
	       xend,
	       r_ns, dr_n, r_ne,
	       z_ns, dz_n, z_ne,
	       phis, dphi, phie,
	       alphas, dalpha, alphae,
	       gammas, dgamma, gammae,
	       FillingTime,
	       CleaningTime,
	       RampUpTime,
	       FullFieldTime,
	       RampDownTime,
	       BFeldSkalGlobal,
	       EFeldSkal);

	// logging parameters to *log.out
	fprintf(LOGSCR, "\nStart parameters: \n"
	                "Energy (min, step, max): %.17LG neV/eV, %.17LG neV/eV, %.17LG neV/eV\n"
	                "Maximum runtime: %.17LG s\n"
	                "r (min, step, max): %.17LG m, %.17LG m, %.17LG m\n"
	                "z (min, step, max): %.17LG m, %.17LG m, %.17LG m\n"
	                "phi (min, step, max): %.17LG m, %.17LG m, %.17LG m\n"
	                "Alpha (min, step, max): %.17LG Grad, %.17LG Grad, %.17LG Grad\n"
	                "Gamma (min, step, max): %.17LG Grad, %.17LG Grad, %.17LG Grad\n"
	                "Filling time %LG s\n"
	                "Cleaning time: %.17LG s\n"
	                "Ramp up time: %.17LG s\n"
	                "Full field time: %.17LG s\n"
	                "RampDownTime: %.17LG s\n"
	                "B field scaling factor: %.17LG s\n"
	                "E field scaling factor: %.17LG s\n",
	                EnergieS, dEnergie, EnergieE,
	                xend,
	                r_ns, dr_n, r_ne,
	                z_ns, dz_n, z_ne,
	                phis, dphi, phie,
	                alphas, dalpha, alphae,
	                gammas, dgamma, gammae,
	                FillingTime,
	                CleaningTime,
	                RampUpTime,
	                FullFieldTime,
	                RampDownTime,
	                BFeldSkalGlobal,
	                EFeldSkal);	
}
//======== end of Startbed ======================================================================================



void ausgabe(long double x2, long double *ystart, long double vend, long double H){
		
	if((x2>=xend)) 
	{                                                           //Zeit abgelaufen
		if (decay)
			kennz=8;   // neutron decayed 
		else 
			kennz=1;   // particle survived until xend			
	}
		
		if(vend>0) gammaend= acosl(ystart[4]/vend) /conv;
		else gammaend=0;
		  alphaend= atan2l(ystart[6]*ystart[1],ystart[2])/conv;
		
	
		// output of end values
		fprintf(ENDLOG,"%i %i %i %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %LG %i %i %LG %LG %li %LG %LG %LG %LG %LG %LG %i %LG %LG %LG\n",
		jobnumber,protneut, polarisation,xstart,r_n,phi_n,z_n,NeutEnergie*1.0e9,v_n,alpha,gammaa,ystart[1],phiend,ystart[3],vend,alphaend,gammaend,x2,H,kennz, NSF,RodFieldMultiplicator, BFflipprob,nrefl, vladmax,vladtotal,thumbmax,
		trajlengthsum,(H-Hstart),Hmax,AbsorberHits, BFeldSkal, EFeldSkal, lossprob);
	
	
		fflush(ENDLOG);
    return;
}


// Diese Funktionen werten die in yp[][] gespeichterten Werte aus und berechnen daraus die
// Neutronenverteilung in der r-z-Ebene. Dies wird �ber die Auswertung der Aufenthalts-
// wahrscheinlichkeit des Neutrons in Quadraten erreicht.
// Wie bei der Interpolation gibt es einen Vektor f�r die r-Werte der Ecken der Quadrate,
// eine f�r die z-Werte und eine Matrix in der die Verteilung abgespeichert wird

// Vorbereitung der r und z-Vektoren
void prepndist(int k)
{
long double dooo, increment;
//v=300;       //  2mm steps in r-direction 0-60cm
//w=1200;      //  z: -0.5 bis1.7
increment = 0.6 / v;  // Schrittweite des Gitters (0.002m)

ndistr=dvector(0,v);       // r-Werte
ndistz=dvector(0,w);       // z-Werte
ndistW=dmatrix(0,v,0,w);   // Verteilung

        // r-Vektor bef�llen, er enth�lt die Ecken der Auswertungsquadrate
        dooo=0.0;
        for(int i=0;i<=v;i++)
        {
				ndistr[i] = dooo;
			  dooo = dooo + increment;
        }

        // z-Vektor bef�llen, er enth�lt die Ecken der Auswertungsquadrate
        dooo=-0.5;
        for(int i=0;i<=w;i++)
        {
				ndistz[i] = dooo;
				dooo = dooo + increment;
        }

        // Matrix der Verteilung auf 0 setzen
        for(int i=0;i<=v;i++)
        {
                for(int j=0;j<=w;j++)
                {
                ndistW[i][j]=0.0;
                }
        }

 }

// Bef�llung der Verteilungsmatrix
void fillndist(int k)
{
	
	int ir1, iz1, ir2, iz2;

	hunt(ndistr, v, yp[1][1], &ir1);     // Index des ersten Punktes suchen
	hunt(ndistz, w, yp[3][1], &iz1);
	//ir1 = (yp[1][1] - conv_rA) / (conv_rB);
	//iz1 = (yp[3][1] - conv_zA) / (conv_zB);
	if((ir1>v)||(iz1>w)) printf("Ndist Error %d %d \n",ir1,iz1);
	
	ir2=ir1; iz2=iz1;

	for(int klauf=2;klauf<=kount;klauf++)
	{
		ir1=ir2; iz1=iz2;
		hunt(ndistr, v, yp[1][klauf], &ir2);   // Index des n�chsten Punktes suchen
		hunt(ndistz, w, yp[3][klauf], &iz2);
		//ir2 = (yp[1][klauf] - conv_rA) / (conv_rB);
		//iz2 = (yp[3][klauf] - conv_zA) / (conv_zB);
//		printf("%i,%i,",ir2,iz2);
		if((ir1>v)||(iz2>w)) printf("Ndist Error %d %d \n",ir2,iz2);
		
		if ((ir2==ir1) && (iz2==iz1))                     // sind die Pkte im gleichen Quadrat?
		{
			ndistW[ir1][iz1]=ndistW[ir1][iz1] + xp[klauf]-xp[klauf-1];            // Zeit zum Quadrat dazuaddieren
		}
		else
		{
	
			ndistW[ir1][iz1]=ndistW[ir1][iz1] + (xp[klauf]-xp[klauf-1])/2;        // H�lfte der Zeit zum ersten Quadrat
			ndistW[ir2][iz2]=ndistW[ir2][iz2] + (xp[klauf]-xp[klauf-1])/2;        // H�lfte der Zeit zum anderen
		}
	}

}



// Ausgabe in ndist.out
void outndist(int k)
{
	//printf("\nOutputting the particle spacial distribution... \n");
		
	ostringstream path;
	path << outpath << "/" << jobnumber << "ndist.out";
	FILE *NDIST=fopen(path.str().c_str(),mode_w);
	int Treffer = 0;
	
	//FILE *NDIST=fopen("ndist.out",mode_w);
	
	fprintf(NDIST,"Rindex Zindex Rmtlpkt Zmtlpkt Whk Treffer\n");
	long double increment = ndistr[1] - ndistr[0];
	int i,j;
	for(i=0;i<=v;i++){
		for(j=0;j<=w;j++){
			if (ndistW[i][j] != 0) Treffer =1;
			if (ndistW[i][j] == 0) Treffer =0;	
			fprintf(NDIST,"%i %i %.5LG %.5LG %.17LG %i\n",i,j,(ndistr[i]+(increment/2.0)),(ndistz[j]+(increment/2.0)), ndistW[i][j], Treffer);
		}
	}
	fclose(NDIST);
	//printf("\nDone outputting the particle spacial distribution! \n");
	
}


void OutputState(long double *y, int l){
	ostringstream stateoutfile;
	stateoutfile << outpath << "/" << jobnumber << "state.out";
	STATEOUT = fopen(stateoutfile.str().c_str(),mode_w);	
	printf("\n \n Something happened! \n \n ");
	fprintf(STATEOUT,"In this file the state of the programm will be written, when something happens for debug porposes... \n");
	fprintf(STATEOUT,"Jobnumber: %i \n", jobnumber);
	fprintf(STATEOUT,"Particle number: %i \n", iMC);
	fprintf(STATEOUT,"Starting energy: %.10LG neV\n", NeutEnergie*1e9);
	fprintf(STATEOUT,"Time: %.10LG s\n", x2);
	fprintf(STATEOUT,"Endtime: %.10LG s\n", xend);
	fprintf(STATEOUT,"Particle ID: %i \n", kennz);
	fprintf(STATEOUT,"Current energy: %.10LG neV\n", H);
	fprintf(STATEOUT,"Spacial coordinates (y[i]): r: %.10LG, phi: %.10LG, z:%.10LG \n", y[1], y[5], y[3]);	
	fprintf(STATEOUT,"Spacial coordinates (ystart[i]): r: %.10LG, phi: %.10LG, z:%.10LG \n", ystart[1], ystart[5], ystart[3]);	
	fprintf(STATEOUT,"Velocities: rdot: %.10LG, phidot: %.10LG, zdot: %.10LG, vabs: %.10LG \n", y[2], y[6], y[4], sqrtl(fabsl(y[2]*y[2]+y[1]*y[1]*y[6]*y[6]+y[4]*y[4])));
	fprintf(STATEOUT,"Magnetic Field: Br: %.10LG, Bphi: %.10LG, Bz: %.10LG, Babs: %.10LG \n", Br,Bphi,Bz,Bws);
	fprintf(STATEOUT,"Magnetic Field Derivatives: \n  %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG %.10LG \n", dBrdr,dBrdphi,dBrdz,dBphidr,dBphidphi,dBphidz,dBzdr,dBzdphi,dBzdz);
	fprintf(STATEOUT,"Reflection state (1 is on): %i \n",reflekt);
	fprintf(STATEOUT,"Diffuse reflection state (1 is specular, 2 diffuse, 3 statistical): %i \n",diffuse);
	fprintf(STATEOUT,"Diffuprob:  %.10LG\n", DiffProb);
	fprintf(STATEOUT,"Rmax: %.10LG\n", rmax);
	fprintf(STATEOUT,"Wanddicke %.10LG\n",wanddicke);
	fprintf(STATEOUT,"Epsi: %.10LG\n", epsi);
	fprintf(STATEOUT,"timestep: %.10LG\n", (x2-x1));
	fprintf(STATEOUT,"Reflekt while Rampdown: %i\n", rdreflekt);
	fprintf(STATEOUT,"Slit: %i\n", slit);
	kennz=88;
	stopall =1;
	exit(-1);
	iMC = MonteCarloAnzahl +1;
}

void csleep(int sec){
time_t start_time, cur_time;

         time(&start_time);
         do
         {
                 time(&cur_time);
         }
         while((cur_time - start_time) < sec);	
}
