/* 	This program generates flamelets solving 
 * 	a 1D counterflow diffusionflame
 *
 *	Author: Federica Ferraro & Hagen Müller
 * 	Universitaet der Bundeswehr Muenchen
 *
 *      email: federica.ferraro@unibw.de

 *      Last Version: 25 April 2013
*/

  #include "cantera/Cantera.h"
  #include "cantera/thermo.h"
  #include "cantera/oneD/Sim1D.h"
  #include "cantera/oneD/Inlet1D.h"
 
  #include "TransportBase_2.h"
  #include "TransportFactory_2.h"
  #include "StFlow_2.h"

  #include "cantera/IdealGasMix.h"
  #include "cantera/equilibrium.h"
  #include "cantera/transport.h"
  #include "cantera/base/plots.h"
  #include <boost/math/special_functions/erf.hpp>
  #include <iostream>
  #include <ostream>
  #include <fstream>
  #include <stdio.h>
  #include <stdlib.h>
  #include <sstream>
  #include <string>
  using namespace Cantera;
  using namespace std;


  void output_TEC(std::ostream& s, const std::string& title,
               const std::vector<std::string>& names,
               const Array2D& data)
  {
    int i,j;

    int npts = static_cast<int>(data.nColumns());
    int nv = static_cast<int>(data.nRows());
    s << "TITLE     = \"" + title + "\"" << endl;
    s << "VARIABLES = " << endl;
    for (i = 0; i < npts; i++) {
        s << " \"" << names[i] << " \" " ;
  }
    s << endl;
    s << "ZONE T=\"zone1\"" << endl;
    s << " I=" << nv << ",J=1,K=1,F=POINT" << endl;
    s << "DT=( ";
    //for (i = 0; i < nv; i++) {
        s << " SINGLE";
    //}
    s << " )" << endl;
    for (i = 0; i < nv; i++) {
        for (j = 0; j < npts; j++) {
            s << data(i,j) << " ";
        }
        s << endl;

    }
    cout << " output TECPLOT File written:" << title << endl;
  }

//****************************************************************
void output_Excel(std::ostream& s, const std::string& title,
                 const std::vector<std::string>& names,
                 const Array2D& data)
{
    int i,j;
    int npts = static_cast<int>(data.nColumns());
    int nv = static_cast<int>(data.nRows());
    for (i = 0; i < npts; i++) {
        s << names[i];
        if (i != npts-1) {
            s << ",";
        }
    }
     s << endl;

     for (i = 0; i < nv; i++) {
        for (j = 0; j < npts; j++) {
            s << data(i,j);
            if (j != npts-1) {
                s << ",";
            }
        }
        s << endl;
   }
cout << " output excel File written:" << title << endl;
}

//****************************************************************************
//****************************************************************************
//*******************************************************************

  void flamelet()
  {
cout<<" "<<endl;
cout<<"**********************************************************************************"<<endl;
cout<<"*                                                                                *"<<endl;
cout<<"*    This program generates flamelets solving 1D counterflow diffusion flame     *"<<endl;
cout<<"*    and it is part of the package flameletFoam.                                 *"<<endl;
cout<<"*                                                                                *"<<endl;
cout<<"*    Author: *Federica Ferraro & Hagen Müller                                    *"<<endl;
cout<<"*    Universitaet der Bundeswehr Muenchen                                        *"<<endl;
cout<<"*                                                                                *"<<endl;
cout<<"*    *email: federica.ferraro@unibw.de                                           *"<<endl;
cout<<"*    *email: hagen.mueller@unibw.de                                              *"<<endl;
cout<<"*                                                                                *"<<endl;
cout<<"*    Last Version: 28 February 2014                                              *"<<endl;
cout<<"*                                                                                *"<<endl;
cout<<"**********************************************************************************"<<endl;  
cout<<" "<<endl; 


//----------Reading the input file--------------------------------------------------------------

    ifstream input_file;
    input_file.open("input.txt", ios::in);
    if (!input_file) {cout<<"Unable to open file input.txt! "<<endl;}

    doublereal pressure;        // Pressure
    doublereal temp_l;          // Temperature flow left side [K]
    doublereal temp_r;          // Temperature flow right side [K]
    doublereal Z_st;            // Stoichiometric mixture fraction
    
    doublereal temp_equil;      // Temperature used in equilibrium calculation
   
    doublereal mdot;
    long int nz;                // Number of grid points
    long int maxnz;             // Maximum number of grid points
    doublereal lz;              // Domain length [m]
    char s[256];
    char composition_l[256];
    char composition_r[256];
    char chem_mech[256];
    char * pch;
    int l=0;
    bool refine_grid;

    while(!input_file.eof())
    {
       l=l+1;
       input_file.getline(s,256);

       if(l==1){
       pressure= atof(s);}
       else if(l==2){
       temp_l=atof(s);}
       else if(l==3){
       temp_r=atof(s);}
       else if(l==4){
       temp_equil=atof(s);}
       else if(l==5){
       mdot=atof(s);}
       else if(l==6){
       Z_st=atof(s);}
       else if (l==7)
       {
	       pch=strchr(s ,' ');
	       strncpy (composition_l, s, pch-s+1 );
	       composition_l[pch-s] = '\0';
       }
       else if (l==8)
       {
	       pch=strchr(s ,' ');
	       strncpy (composition_r, s, pch-s+1 );
	       composition_r[pch-s] = '\0';
       }
       else if (l==9)
       {
    	    pch=strchr(s ,' ');
    	    strncpy (chem_mech, s, pch-s+1 );
    	    chem_mech[pch-s] = '\0';
       }
       else if (l==10){
       lz=atof(s);}
       else if (l==11){
       nz=atof(s);}
       else if (l==12){
       maxnz=atoi(s);
       input_file>>refine_grid;}
    }

    input_file.close();

    doublereal mdot_l=mdot;          // Mass flow left side (density * velocity) [kg/m^2/s]
    doublereal mdot_r=mdot;          // Mass flow right side (density * velocity)[kg/m^2/s]

    int loglevel=1;
    unsigned int i, h,k;

    IdealGasMix gas(chem_mech); // Write in brackets the path of the kinetics mechanism file
    unsigned int nsp=gas.nSpecies();
        
    vector_fp x;
    x.resize(nsp);

    int flowdomain=1;
    
    //----------------------------------------------------------------------------
    //----------Create the flow -----------------------------------------------
    AxiStagnFlow flow(&gas);

    //--------- Create an initial grid-----------------------------------------
    doublereal* z=new double[nz+1];
    doublereal  dz=lz/((doublereal)(nz-1));
    int iz;                          // Index for the grid points  
    
    for ( iz=0; iz<nz; iz++) {
          z[iz]=((doublereal)iz)*dz;  
        }
 
    flow.setupGrid(nz, z);

    //-------Choose the Transport Model:
    Transport* tr_lew = newTransportMgr("Lewis1", &gas);
    flow.setTransport(*tr_lew);
    flow.setKinetics(gas);
     
    // --------Set the state of the two inlets---------------------------------

    //---------Left side-------------------------------------------------------
    Inlet1D left_side;

     //----Set the initial composition in mole fractions
    left_side.setMoleFractions(composition_l);

    //-----y_l and y_r are the vectors with initial boundary conditions, used to set an initial solution
    double* y_l=new double[nsp];
    
    char curSpecie[10];
    char curMoleFr[10];
    char * pch1;
    char * pch2;
    char * pch3;
    pch1=strchr(composition_l ,':');
    pch2=strchr(composition_l ,',');
    pch3=strchr(composition_l ,'\0');

    memmove (curSpecie, composition_l, pch1-composition_l+1);
    curSpecie[pch1-composition_l] = '\0';

    while (pch1!=NULL)
    {
        if (pch2!=NULL)
        {
        	memmove (curMoleFr, composition_l+(pch1-composition_l)+1, (pch2-composition_l)-(pch1-composition_l)-1);
            curMoleFr[(pch2-composition_l)-(pch1-composition_l)-1] = '\0';
            for ( h=0; h<nsp; h++)
            {
      	         if (gas.speciesIndex(curSpecie)==h) {y_l[h]=atof(curMoleFr);}
            }
            pch1=strchr(pch1+1,':');
            memmove (curSpecie, composition_l+(pch2-composition_l)+1, (pch1-composition_l)-(pch2-composition_l)-1);
            curSpecie[(pch1-composition_l)-(pch2-composition_l)-1] = '\0';
            pch2=strchr(pch2+1,',');
        }
        else
        {
        	memmove (curMoleFr, composition_l+(pch1-composition_l)+1, (pch3-composition_l)-(pch1-composition_l)-1);
            curMoleFr[(pch3-composition_l)-(pch1-composition_l)-1] = '\0';
            for ( h=0; h<nsp; h++)
            {
      	         if (gas.speciesIndex(curSpecie)==h) {y_l[h]=atof(curMoleFr);}
            }
            pch1=strchr(pch1+1,':');
        }
    }
     
    left_side.setMdot(mdot_l);
    left_side.setTemperature(temp_l);

    //--------right side--------------------------------------------------------
    Inlet1D right_side;

    //----Set the initial composition in mole fractions
    right_side.setMoleFractions(composition_r);

    double* y_r=new double[nsp];

    pch1=strchr(composition_r ,':');
    pch2=strchr(composition_r ,',');
    pch3=strchr(composition_r ,'\0');

    memmove (curSpecie, composition_r, pch1-composition_r+1);
    curSpecie[pch1-composition_r] = '\0';

    while (pch1!=NULL)
    {
        if (pch2!=NULL)
        {
        	memmove (curMoleFr, composition_r+(pch1-composition_r)+1, (pch2-composition_r)-(pch1-composition_r)-1);
            curMoleFr[(pch2-composition_r)-(pch1-composition_r)-1] = '\0';
            for ( h=0; h<nsp; h++)
            {
      	         if (gas.speciesIndex(curSpecie)==h){y_r[h]=atof(curMoleFr);}
            }
            pch1=strchr(pch1+1,':');
            memmove (curSpecie, composition_r+(pch2-composition_r)+1, (pch1-composition_r)-(pch2-composition_r)-1);
            curSpecie[(pch1-composition_r)-(pch2-composition_r)-1] = '\0';
            pch2=strchr(pch2+1,',');
        }
        else
        {
        	memmove (curMoleFr, composition_r+(pch1-composition_r)+1, (pch3-composition_r)-(pch1-composition_r)-1);
            curMoleFr[(pch3-composition_r)-(pch1-composition_r)-1] = '\0';
            for ( h=0; h<nsp; h++)
            {
      	         if (gas.speciesIndex(curSpecie)==h) {y_r[h]=atof(curMoleFr);}
            }
            pch1=strchr(pch1+1,':');
        }
    }

    right_side.setMdot(mdot_r);
    right_side.setTemperature(temp_r);
  

    // --------domain definition
     
    std::vector<Domain1D*> domains;
    domains.push_back(&left_side);
    domains.push_back(&flow);
    domains.push_back(&right_side);   

    Sim1D flame(domains);

    cout<< "Boundary conditions: Left side"<<endl;
    left_side.showSolution(0);
    cout<< "Boundary conditions: Right side"<<endl;
    right_side.showSolution(0);
    
    int sol_ini;
    doublereal new_lz, scale_factor;

    cout << "Do you want to use the initial solution from file? Type 1 for yes or 0 to no: " << endl;
    cin >> sol_ini;
    
    if (sol_ini == 1)
    { 
       flame.restore("solution.xml","id_01");

       nz=flow.nPoints();
       doublereal* newz=new double[nz+1];
       if (flow.grid(nz-1)!=lz)
       {
           cout<<" "<<endl;
           cout<<"*****************************************************************************************"<<endl;
           cout<< endl << "The domain dimension in the solution is different from the value in the input file";
           cout<< endl << "The domain in the solution file is "<< flow.grid(nz-1)<<" m";
           cout<< endl << "The initial solution is scaled to the new domain size: "<< lz <<" m"<<endl;
           cout<<"*****************************************************************************************"<<endl;
           cout <<" " <<endl;

           new_lz=lz;            //update domain dimension
           scale_factor = new_lz/flow.grid(nz-1);
           for ( iz=0; iz<nz; iz++)
           {
                newz[iz] = flow.grid(iz) * scale_factor;
           }
           flow.setupGrid(nz, newz);
       }
    }

    else
    {
        // In this part the equilibrium solution is provided as an initial guess
        double phi = 0.0;
        string fuel;
        string oxidizer;
        cout<<"*********************************************************************************"<<endl;
        cout << "An equilibrium solution will be used as initial guess"<< endl;
        cout << "Please enter a value for the equivalence ratio phi (0 - 1): "<<endl;
        cin >> phi;
        cout << "Please enter the name of the Fuel (CH4, H2): "<<endl;
        cin >> fuel;
        cout << "Please enter the name of the Oxidizer (O2, air): "<<endl;
        cin >> oxidizer;
        cout<<"*********************************************************************************"<<endl;

        if (fuel == "CH4" && oxidizer == "air")
        {
        	doublereal fa_stoic = 0.105042;
            for ( h=0; h<nsp; h++)
            {
               if (h==gas.speciesIndex("CH4"))
               {
                  x[h]=1.0;
               }
               else if (h==gas.speciesIndex("O2"))
               {
                  x[h]=0.21/phi/fa_stoic;
               }
               else if (h==gas.speciesIndex("N2"))
               {
                  x[h]=0.79/phi/fa_stoic;
               }
               else
               {
                  x[h]=0.0;
               }
            }
        }
        else if (fuel == "CH4" && oxidizer == "O2")
        {
        	doublereal fa_stoic = 0.5;
            for ( h=0; h<nsp; h++)
            {
               if (h==gas.speciesIndex("CH4"))
               {
                  x[h]=1.0;
               }
               else if (h==gas.speciesIndex("O2"))
               {
                  x[h]=1.0/phi/fa_stoic;
               }
               else
               {
                  x[h]=0.0;
               }
            }
        }
        else if (fuel == "H2" && oxidizer == "O2")
        {
        	doublereal fa_stoic = 2;
            for ( h=0; h<nsp; h++)
            {
               if (h==gas.speciesIndex("h2"))
               {
                  x[h]=1.0;
               }
               else if (h==gas.speciesIndex("o2"))
               {
                  x[h]=1.0/phi/fa_stoic;
               }
               else
               {
                  x[h]=0.0;
               }
            }
        }
        else if (fuel == "H2" && oxidizer == "air")
        {
        	doublereal fa_stoic = 0.3471;
            for ( h=0; h<nsp; h++)
            {
               if (h==gas.speciesIndex("H2"))
               {
                  x[h]=1.0;
               }
               else if (h==gas.speciesIndex("O2"))
               {
                  x[h]=0.21/phi/fa_stoic;
               }
               else if (h==gas.speciesIndex("N2"))
               {
                  x[h]=0.79/phi/fa_stoic;
               }
               else
               {
                  x[h]=0.0;
               }
            }
        }
        else
        {
            cout<<"*********************************************************************************"<<endl;
        	cout << "Your fuel/oxidizer selection is not available"<< endl;
        	cout << "Please add your combination in the source code"<< endl;
        	cout<<"*********************************************************************************"<<endl;
        }

        gas.setState_TPX(temp_equil, pressure, DATA_PTR(x));
        double* yin=new double[nsp];
        gas.getMoleFractions(yin);
      
        //----Equilibrium  composition calculation for constant enthalpy and pressure (HP)----------
        try
        {
          equilibrate(gas,"HP");
        }
        catch (CanteraError& err)
        {
          std::cout << err.what() << std::endl;
        }
       
        double* yout=new double[nsp];
        gas.getMoleFractions(yout);
        doublereal Tad=gas.temperature();
        cout<<"Successful Equilibrium Calculation"<<endl;
        cout <<"phi  "<< phi<<' '<<"Adiabatic Temperature "<< Tad << " K" <<endl;

        //----Set initial profiles of temperature and species from the equilibrium calculation-----------
        vector_fp locs;
        vector_fp value;
        locs.resize(3);
        value.resize(3);
        locs[0]=0;
        locs[1]=0;
        locs[2]=1.0;

        value[0]=temp_l;
        value[1]=Tad;
        value[2]=temp_r;
        flame.setInitialGuess("T",locs,value);

        for (i=0; i<nsp;i++)
        {
           value[0]=y_l[i];
           value[1]=yout[i];
           value[2]=y_r[i];
           flame.setInitialGuess(gas.speciesName(i), locs, value);
        }
    }      
   
    //-----------------------Set parameters for grid refinement-----------------------------------
    double ratio=200.0;
    double slope=0.1;
    double curve=0.2;
    double prune=0.02;

    int log_num;
 
    if(refine_grid==true)
    {
        flame.setMaxGridPoints(flowdomain, maxnz);
        cout<<"*********************************************************************************"<<endl;
        cout<<"Grid Refinement Activated!"<<endl;
        cout<<"*********************************************************************************"<<endl;
    }

    // Start the Calculation
    flow.setPressure(pressure);  
    flow.solveEnergyEqn(); 
    flame.setRefineCriteria(flowdomain,ratio,slope,curve,prune);
    cout<<"Flamelet Calculation Started"<<endl;
    flame.solve(loglevel, refine_grid);
   
    int np=flow.nPoints();
    Array2D solution(np,nsp+6,0.0);
    vector_fp grid;
    
    if (np!= nz) {
    grid.resize(np);}
    else {
    grid.resize(nz);
    }
   
    for ( iz=0; iz<np; iz++) {
         grid[iz]=flow.grid(iz);
    }
    
    // Save the solution
    int n;  
    for (n=0; n<np; n++)
    {
         solution(n,0)=grid[n];
    }
 

    double argonNorm = max(flame.value(flowdomain, flow.componentIndex("AR"),0),flame.value(flowdomain, flow.componentIndex("AR"),np-1));
    for (n=0; n<np; n++)
    {
       for (i=1; i<5;i++)
       {
          solution(n,i)=flame.value(flowdomain, i-1,n);
       }

       for (i=5; i<nsp+5;i++)
       {
          if (flame.value(flowdomain, i-1,n) > 0.0){solution(n,i) = flame.value(flowdomain, i-1,n);}
          else if (flame.value(flowdomain, i-1,n) < 0.0)
          {
              solution(n,i) = 0.0;
              flame.setValue(flowdomain, i-1,n,0.0);
          }
       }

       // Calculate the mixture fraction from Argon
       solution(n,nsp+5)=flame.value(flowdomain, flow.componentIndex("AR"),n)/argonNorm;
    }

    std::vector<std::string> names;
    names.resize(nsp+6);
    names[0]= "x"; 

    names[1]  = "u";
    names[2]  = "V";
    names[3]  = "T";
    names[4]  = "Lambda";
    for ( k=0; k<nsp; k++){
    names[5+k]=gas.speciesName(k);
    }
    names[nsp+5]  = "Z";

    //-------Scalar Dissipation Rate calculation
    double a = sqrt(abs(solution(1,4))/ gas.density());
    double dissRate = 2.0*a/M_PI*exp(-2.0*pow(boost::math::erfc_inv(2.0*Z_st),2));

    cout << "Scalar Dissipation Rate: " << dissRate << " 1/s"<< endl;

    //-------Construct Filename
    ostringstream strs;
    strs << dissRate;
    std::string str_dissRate = strs.str();
    string strfile = "canteraTables/canteraTable_" + str_dissRate + ".csv";

    char * filename = new char[strfile.length()];
    strcpy(filename,strfile.c_str());

    // Print the flamelet table to be used in OPENFOAM
    ofstream myfile;
    myfile.open(filename);
    output_Excel(myfile, filename, names, solution);
 
    string title;
    string id;
    string desc;
    cout << "Do you want to write the initial guess file? press 1 for yes or 0 to no";
    cin >> log_num;

    if(log_num==1){
    remove ("solution.xml"); 
    flame.save(title= "solution.xml",id="id_01", desc="case_1");}

  }

  int main()
  {
    try
    {
       flamelet();
    }
    catch (CanteraError& err)
    {
       std::cout << err.what() << std::endl;
    }
  }
