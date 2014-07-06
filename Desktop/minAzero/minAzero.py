#!/usr/bin/python

#####################################################
#
#  10 January 2011
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#  This script builds C++ code to minimize a path integral
#   formulation of a parameter estimation problem with
#   the optimization software IPOPT
#
#  Specifically, given a vector-field (model) of the form:
#
#           dx_1(t) = G_1(x_1(t),x_p(t),q)
# 
#           dx_p(t) = G_p(x_1(t),x_p(t),q)
#
#          where x_p denotes 1 or more equations in the model,
#
#  this code takes the discAzerod vector field and objective 
#  function (discAzerod in companion script discAzero.py),
#  and builds the requisite IPOPT functions to solve the 
#  resulting optimization problem.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but could
#  easily be modified for use with other optimization software.
#
######################################################


#  For ease of use, all necessary C++ files are built with one keyboard
#   command.  The following four scripts each write a necessary C++ file.
#   See the individual scripts for more information.


import discAzero

import makecppAzero
import makehppAzero
import makemakeAzero
import makeoptAzero
# DiscAzero.py reads the file equations.txt and sets up the given vector field
#  in the correct format.
print "Writing the problem files\n"
# Import the problem name and change to upper and lower case
prob = discAzero.Problem
probu = prob.upper()
probl = prob.lower()

FILE = probl + 'minAzero_nlp.cpp'

nU = discAzero.nU
nP = discAzero.nP
nY = discAzero.nY
nI = discAzero.nI
nF = discAzero.nF
nM = discAzero.nM
# The name of the IPOPT file to be written to
f = open(FILE,'w')

# The following write commands are writing C++ code.

# Front matter
f.write('// %s.cpp\n' % probl)
f.write('// Nonlinear Ipopt program\n')
f.write('\n// Author: Bryan A. Toth\n// btoth@physics.ucsd.edu\n\n')
f.write('#include \"%sminAzero_nlp.hpp\"\n' % probl)
f.write('#include <cmath>\n\
#include <cstdio>\n\
#include <iostream>\n\
#include <fstream>\n\
#include <string>\n\
#include <stdlib.h>\n\
#include <cstring>\n')

f.write('#ifdef HAVE_CSTDIO\n\
# include <cstdio>\n\
# include <iostream>\n\
# include <fstream>\n\
# include <string>\n\
# include <stdlib.h>\n\
#else\n\
# ifdef HAVE_STDIO_H\n\
#  include <stdio.h>\n\
# else\n\
#  error \"don\'t have header file for stdio\"\n\
# endif\n\
#endif\n\n')

f.write('using namespace Ipopt;\n\n')
f.write('using namespace std;\n\n')
for i in range(nF):
  args = discAzero.Funcarg[i]
  f.write('double %s(' % discAzero.Sf[i])
  for j in range(args-1):
     f.write('double, ')
  f.write('double);\n')
  f.write('double %sjac(' % discAzero.Sf[i])
  for j in range(args-1):
     f.write('double, ')
  f.write('double, int);\n')
  f.write('double %shes(' % discAzero.Sf[i])
  for j in range(args-1):
     f.write('double, ')
  f.write('double, int, int);\n')
f.write('// constructor\n\
%s_NLP::%s_NLP()\n\
{\n' % (probu, probu))

# Define problem parameters
f.write('  nU=%d;\n' % nU)
f.write('  nP=%d;\n' % nP)
f.write('  nY=%d;\n' % nY)
f.write('  nM=%d;\n' % nM)
f.write('  nI=%d;\n\n' % nI)

# Define variable names 
f.write('\
     K11val = new double[nU];\n\
     K11val2 = new double[nU];\n\
     K11valp1 = new double[nU];\n\
     dK11val = new double[nU];\n\
     dK11val2 = new double[nU];\n\
     dK11valp1 = new double[nU];\n\
     Xdval = new double[nM];\n\
     Xdval2 = new double[nM];\n\
     Xdvalp1 = new double[nM];\n')
f.write('\
     Xval = new double[nY];\n\
     Xval2 = new double[nY];\n\
     Xvalp1 = new double[nY];\n')
f.write('\
     Pval = new double[nP];\n')
f.write('\
     Ival = new double[nI];\n\
     Ival2 = new double[nI];\n\
     Ivalp1 = new double[nI];\n')
# Commands to read specs.txt file
f.write('\
     string buffer;\n\
     specs = new string[6+nP+nY+nI+2*nU+nM];\n\
     \n\
     int count;\n\
     count = 0;\n\
     \n\
     ifstream fin ("specs.txt");\n')
f.write("     if (fin.is_open())\n\
     {\n\
       while (! fin.eof())\n\
       {\n\
         getline (fin,buffer);\n\
	 if (buffer[0] !='#')\n\
	   {\n\
	   specs[count] = buffer;\n\
	   count++;\n\
	   }\n\
       }\n\
       fin.close();\n\
     }\n\
     \n\
     else cout << \"Unable to open file\";\n")

# Write Time to file
#  Time is a misnomer - this is a measure of the the number
#  of time steps used in the problem
f.write('    Time = atoi(specs[0].c_str());\n')
# Write skip to file
#  Skip is a dummy variable to allow the use of various parts
#  of a given data file
f.write('    skip = atoi(specs[1].c_str());\n')
# Write hstep to file
#  Hstep is the time-step of the discretization
f.write('    hstep = atof(specs[2].c_str());\n\n')

# Write open data file to file
f.write('    string filename;\n')
f.write('    int ret;\n')

# Data for each variable that is being coupled in to the vector field
for i in range(nM):
    f.write('\
    %s = new double[2*Time+1];\n\
    %s = new double[skip];\n\n\
    FILE *pFile%d;\n' % (discAzero.Ldata[i],discAzero.Ldata[i]+'dummy',i))
    f.write('\
    filename = specs[%d];\n' % (3+i))
    f.write('\
    pFile%d = fopen(filename.c_str(),"r");\n' % i)

# Read data from data file
    temp1 = "%lf"
    f.write('\n\
    for(Index jt=0;jt<skip;jt++)\n\
	{\n\
	ret = fscanf (pFile%d, "%s", &%s[jt]);\n\
	if (ret == EOF) break;\n\
	}\n\
    for(Index jt=0;jt<2*Time+1;jt++)\n\
	{\n\
	ret = fscanf (pFile%d, "%s", &%s[jt]);\n\
	if (ret == EOF) break;\n\
	}\n\
    fclose (pFile%d);\n' % (i,temp1, discAzero.Ldata[i]+'dummy',i,temp1,discAzero.Ldata[i],i))
#########  End for loop #############

# Open data file for stimulus
for i in range(nI):
    f.write('\
    %s = new double[2*Time+1];\n\
    %s = new double[skip];\n\n\
    FILE *qFile%d;\n' % (discAzero.Lstimuli[i],discAzero.Lstimuli[i]+'dummy',i))
    f.write('\
    filename = specs[%d];\n' % (3+nM+i))
    f.write('\
    qFile%d = fopen(filename.c_str(),"r");\n' % i)
# Read stimuli data
    temp1 = "%lf"
    f.write('\n\
    for(Index jt=0;jt<skip;jt++)\n\
        {\n\
	ret = fscanf (qFile%d, "%s", &%s[jt]);\n\
	if (ret == EOF) break;\n\
        }\n\
    for(Index jt=0;jt<2*Time+1;jt++)\n\
        {\n\
	ret = fscanf (qFile%d, "%s", &%s[jt]);\n\
	if (ret == EOF) break;\n\
	}\n\
    fclose (qFile%d);\n' % (i,temp1,discAzero.Lstimuli[i]+'dummy',i,temp1,discAzero.Lstimuli[i],i))
#########  End for loop #############

# Read in the initial and boundary conditions for all variables into arrays
f.write('\
    int rows = nY+2*nU+nP;\n\
    bounds = new double*[rows];\n\
    for (Index i=0;i<rows;i++) bounds[i] = new double[4];\n\
    int toggle=0;\n\
    if (specs[3+nM+nI] == "1") toggle = 1;\n\
    int counter;\n\
    for(Index k=0;k<rows;k++)\n\
       {\n\
       counter=0;\n\
       char* tmp = new char[specs[4+nM+nI+toggle+k].size()+1];\n\
       strcpy( tmp, specs[4+nM+nI+toggle+k].c_str() );\n\
       char *ptr = strtok(tmp,",");\n\
       bounds[k][3] = 0.0;\n\
       while(ptr != 0) {\n\
          if(counter<3) {\n\
	     bounds[k][counter] = atof(ptr);\n\
	     }\n\
          if(counter==3) {\n\
             bounds[k][counter] = atof(ptr);\n\
             }\n\
	  ptr = strtok(0,",");\n\
	  counter++;\n\
          }\n\
    }\n\n')
# If initial conditions are in a data file, read in the data file
f.write('\
    if (specs[3+nM+nI] == "1")\n\
       {\n\
       filename = specs[4+nM+nI];\n\
       }\n\n')
  

f.write('\
}\n\n')

f.write('// destructor\n\
%s_NLP::~%s_NLP()\n\
{\n\
  delete [] K11val;\n\
  delete [] K11val2;\n\
  delete [] K11valp1;\n\
  delete [] dK11val;\n\
  delete [] dK11val2;\n\
  delete [] dK11valp1;\n\
  delete [] Xdval;\n\
  delete [] Xdval2;\n\
  delete [] Xdvalp1;\n\
  delete [] Xval;\n\
  delete [] Xval2;\n\
  delete [] Xvalp1;\n\
  delete [] Pval;\n\
  delete [] Ival;\n\
  delete [] Ival2;\n\
  delete [] Ivalp1;\n\
  delete [] specs;\n' % (probu, probu))
for i in range(nM):
    f.write('\
    delete [] %s;\n\
    delete [] %s;\n' % (discAzero.Ldata[i],discAzero.Ldata[i]+'dummy'))
for i in range(nI):
    f.write('\
    delete [] %s;\n\
    delete [] %s;\n' % (discAzero.Lstimuli[i],discAzero.Lstimuli[i]+'dummy'))
f.write('\
  int rows = nY+2*nU+nP;\n\
  for (Index i=0;i<rows;i++) delete [] bounds[i];\n\
  delete [] bounds;\n')
f.write('\n\
}\n\n')
# Start to write individual functions

# GET_NLP_INFO


f.write('// returns the size of the problem\n\
bool %s_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,\n\
				Index& nnz_h_lag, IndexStyleEnum& index_style)\n\n' % probu)
# Number of variables
alpha = 2*(nY+2*nU)
beta = nY+2*nU+nP
f.write('{\n\
  // Number of variables\n\
  n = %d*Time+%d;\n' % (alpha,beta))

# Number of equality constraints
f.write('\n\
  // Number of equality constraints\n\
  m = 0;\n')

# Number of Jacobian nonzero entries
f.write('\n\
  // Number of Jacobian nonzero entries\n\
  nnz_jac_g = 0;\n')

# Number of Hessian non-zeros in lower left of diagonal
# Determine number of large, medium, small elements in Hessian
Hessize = discAzero.Hessize
Hesindex = discAzero.Hesindex
cat = len(discAzero.Smod)
sma = 0
med = 0
lar = 0
for i in range(cat):
  for j in range(i+1):
      if Hesindex[i][j] == -1:
         temp = Hessize[i][j]
         if temp == 's':
            sma = sma + 1
         elif temp == 'm':
            med = med + 1
         elif temp == 'l':
            lar = lar + 1 

f.write('\n\
  // Number of Hessian nonzero entries\n\
  nnz_h_lag = %d*(Time+1)+%d*Time+%d;\n' % (lar, med, sma))

f.write('\n\
  // use the C style indexing (0-based)\n\
  index_style = TNLP::C_STYLE;\n\n\
  return true;\n\
}\n\n\n')



# GET_BOUNDS_INFO


f.write('// returns the variable bounds\n\
bool %s_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,\n\
                                Index m, Number* g_l, Number* g_u)\n\
{\n\
  // Here, the n and m we gave IPOPT in get_nlp_info are passed back to us.\n\
  // If desired, we could assert to make sure they are what we think they are.\n' % probu)
f.write('  assert(n == %d*Time+%d);\n' % (alpha, beta))
f.write('  assert(m == 0);\n\n')

# This takes the bounds information read in from specs.txt, and puts
#  it into the correct spot in a bounds array for all time points.

f.write('  for(Index jt=0;jt<Time+1;jt++) {\n')
f.write('     for(Index var=0;var<nY;var++) {\n')
f.write('        // Bounds for x\n')
f.write('        x_l[(Time+1)*var+jt]=bounds[var][0];\n')
f.write('        x_u[(Time+1)*var+jt]=bounds[var][1];\n')
f.write('        // Bounds for midpoints\n')
f.write('        if(jt<Time) {\n\
       x_l[(Time+1)*(nY+2*nU)+Time*var+jt]=bounds[var][0];\n\
       x_u[(Time+1)*(nY+2*nU)+Time*var+jt]=bounds[var][1];\n\
       }\n\
    }\n')
  #this is where constraint bounds are set
f.write('     for(Index con=0;con<2*nU;con++) {\n')    
f.write('       // Bounds for k\n')
f.write('       x_l[(Time+1)*(nY+con)+jt]=bounds[nY+con][0];\n')
f.write('       x_u[(Time+1)*(nY+con)+jt]=bounds[nY+con][1];\n')
f.write('       // Bounds for midpoints\n')
f.write('       if(jt<Time) {\n\
          x_l[(Time+1)*(nY+2*nU)+Time*(nY+con)+jt]=bounds[nY+con][0];\n\
	  x_u[(Time+1)*(nY+2*nU)+Time*(nY+con)+jt]=bounds[nY+con][1];\n\
	  }\n\
     }\n\n')
f.write('  } // End for loop\n\n')
f.write('     for(Index par=0;par<nP;par++) {\n')
f.write('        // Bounds for parameters\n')
f.write('        x_l[2*Time*(nY+2*nU)+nY+2*nU+par]=bounds[nY+2*nU+par][0];\n')
f.write('        x_u[2*Time*(nY+2*nU)+nY+2*nU+par]=bounds[nY+2*nU+par][1];\n\
              }\n\n')

f.write('  return true;\n\
}\n\n\n' )


# GET_STARTING_POINT

f.write('// returns the initial point for the problem\n\
bool %s_NLP::get_starting_point(Index n, bool init_x, Number* x,\n\
                                   bool init_z, Number* z_L, Number* z_U,\n\
                                   Index m, bool init_lambda,\n\
                                   Number* lambda)\n\
{\n\
  assert(init_x == true);\n\
  assert(init_z == false);\n\
  assert(init_lambda == false);\n\n\
  for (Index i=0; i<n; i++) {\n\
        x[i] = 0.0;\n\
      }\n\n' % probu)
f.write('\
      int skipROWS = skip;\n\
      int ROWS = 2*Time+1;\n\
      int COLS = nY;\n\
    \n\
      double **skipinit = new double* [skipROWS];\n\
      double **init = new double* [ROWS];\n\
      for(Index i=0;i<skipROWS;i++) skipinit[i] = new double[COLS];\n\
      for(Index i=0;i<ROWS;i++) init[i] = new double[COLS];\n\
      \n\
      string filename;\n\
      filename = specs[4+nM+nI];\n\n')

# To start from an initial guess given in a separate data file:
#  Read in the data file
temp1 = "%lf"
f.write('\
      if (specs[3+nM+nI] =="1")\n\
      {\n\
      FILE *initFILE;\n\
      int ret;\n\
      initFILE = fopen(filename.c_str(),"r");\n\
    \n\
      for(Index jt=0;jt<skip;jt++)\n\
          {\n\
	  ret = fscanf (initFILE,"')
for i in range(nY):
   f.write('%s ' % temp1)
f.write('"')
for i in range(nY):
   f.write(',&skipinit[jt][%d]' % i)
f.write(');\n\
	  if (ret == EOF) break;\n\
          }\n\
      for(Index jt=0;jt<2*Time+1;jt++)\n\
          {\n\
	  ret = fscanf (initFILE,"')
for i in range(nY):
   f.write('%s ' % temp1)
f.write('"')
for i in range(nY):
   f.write(',&init[jt][%d]' % i)
f.write(');\n\
	  if (ret == EOF) break;\n\
          }\n\
    fclose (initFILE);\n\
    }\n\n')

# Set the initial starting point into the x[] array, either from
#  the numbers given in specs.txt or from an initial data file

f.write('  for(Index jt=0;jt<Time+1;jt++) {\n')
f.write('     for(Index var=0;var<nY;var++) {\n')    
f.write('       // Initial conditions for x\n')
f.write('       if (specs[3+nM+nI] == "1")\n\
                  {\n\
		  x[(Time+1)*var+jt] = init[2*jt][var];\n\
		  }\n\
	        else\n\
	          {\n\
		  x[(Time+1)*var+jt] = bounds[var][2];\n\
		  }\n')
f.write('       // Initial conditions for midpoints\n')
f.write('       if(jt<Time) {\n\
                  if (specs[3+nM+nI] == "1")\n\
		    {\n\
		    x[(Time+1)*(nY+2*nU)+Time*var+jt] = init[2*jt+1][var];\n\
		    }\n\
		  else\n\
		    {\n\
		    x[(Time+1)*(nY+2*nU)+Time*var+jt] = bounds[var][2];\n\
		    }\n\
		  }\n\
		}\n')

f.write('     for(Index cup=0;cup<2*nU;cup++) {\n')
f.write('       // Initial conditions for k\n')
f.write('       x[(Time+1)*(cup+nY)+jt]=bounds[cup+nY][2];\n')
f.write('       // Initial conditions for midpoints\n')
f.write('       if(jt<Time) {\n\
                  x[(Time+1)*(nY+2*nU)+Time*(cup+nY)+jt]=bounds[cup+nY][2];\n\
	          }\n\
	      }\n')

f.write('  } // End for loop\n\n')

f.write('     for(Index par=0;par<nP;par++) {\n')
f.write('     // Initial conditions for p%d\n' % (i+1))
f.write('     x[2*Time*(nY+2*nU)+nY+2*nU+par]=bounds[nY+2*nU+par][2];\n\
              }\n\n')
f.write('  for(Index i=0;i<ROWS;i++) delete [] init[i];\n\
   delete [] init;\n') 
f.write('  return true;\n\
}\n\n\n') 



# EVAL_F
# Subroutine to calculate the objective value
# Here, and in the following subroutines, strings from the result
#  of symbolic discretization and differentiation in discAzero.py
#  are inserted to the code.

f.write('// returns the value of the objective function\n\
bool %s_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)\n\
{\n' % probu)
f.write('  assert(n == %d*Time+%d);\n' % (alpha, beta))
f.write('  obj_value = 0;\n\n')
f.write('  for(Index jt=0;jt<Time;jt++) {\n\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[jt + i*(Time+1)];\n')
f.write('        Xvalp1[i] = x[jt + i*(Time+1) + 1];\n')
f.write('        Xval2[i] = x[(Time+1)*(nY+2*nU) + jt + i*(Time)];\n')
f.write('     } //end for loop\n\n')

f.write('\n')
f.write('     for(Index i=0;i<nU;i++) {\n')
f.write('        K11val[i] = x[jt + nY*(Time+1) + 2*i*(Time+1)];\n')
f.write('        K11valp1[i] = x[jt + nY*(Time+1) + 2*i*(Time+1) + 1];\n')
f.write('        K11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i)*Time + jt];\n')
f.write('        dK11val[i] = x[jt + (nY+2*i+1)*(Time+1)];\n')
f.write('        dK11valp1[i] = x[jt + (nY+2*i+1)*(Time+1)+1];\n')
f.write('        dK11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i+1)*Time + jt];\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*jt];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdval2[%d] = %s[2*jt+1];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdvalp1[%d] = %s[2*jt+2];\n' % (i,discAzero.Ldata[i]))

for i in range(nI):
    f.write('     Ival[%d] = %s[2*jt];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ival2[%d] = %s[2*jt+1];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ivalp1[%d] = %s[2*jt+2];\n' % (i,discAzero.Lstimuli[i]))

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')
f.write('    obj_value += %s; \n\n' % (discAzero.strObj[0]))
for i in range(nY):
   f.write('    obj_value += bounds[%s][3]*(%s); \n\n' % (i,discAzero.strObj[1][i]))
  
f.write('  } //end for loop\n\n')

# Add code for last element

f.write('// Add last element\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[Time + i*(Time+1)];\n')
f.write('        Xvalp1[i] = 0;\n')
f.write('        Xval2[i] = 0;\n')
f.write('     } //end for loop\n\n')

f.write('\n')
f.write('     for(Index i=0;i<nU;i++) {\n')
f.write('        K11val[i] = x[Time + nY*(Time+1) + 2*i*(Time+1)];\n')
f.write('        K11valp1[i] = 0;\n')
f.write('        K11val2[i] = 0;\n')
f.write('        dK11val[i] = x[Time + (nY+2*i+1)*(Time+1)];\n')
f.write('        dK11valp1[i] = 0;\n')
f.write('        dK11val2[i] = 0;\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*Time];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdval2[%d] = 0;\n' % i)
    f.write('     Xdvalp1[%d] = 0;\n' % i)

for i in range(nI):
    f.write('     Ival[%d] = %s[2*Time];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ival2[%d] = 0;\n' % i)
    f.write('     Ivalp1[%d] = 0;\n' % i)

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')

f.write('  obj_value += %s;\n\n' % (discAzero.strObj[0]))

# Adding a line to divide the overall objective function by 2T + 1
# This normalizes the objective function

f.write('  obj_value = obj_value/(2*Time+1);\n\n')

f.write('  return true;\n\
}\n\n\n')



# EVAL_GRAD_F

f.write('// return the gradient of the objective function grad_{x} f(x)\n\
bool %s_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)\n\
{\n' % probu)
f.write('  assert(n == %d*Time+%d);\n\n' % (alpha, beta))

f.write('  for(Index i=0;i<n;i++) {\n')
f.write('     grad_f[i] = 0;\n')
f.write('  }\n\n')

f.write('  for(Index jt=0;jt<Time;jt++) {\n\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[jt + i*(Time+1)];\n')
f.write('        Xvalp1[i] = x[jt + i*(Time+1) + 1];\n')
f.write('        Xval2[i] = x[(Time+1)*(nY+2*nU) + jt + i*(Time)];\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index i=0;i<nU;i++) {\n')
f.write('        K11val[i] = x[jt + nY*(Time+1) + 2*i*(Time+1)];\n')
f.write('        K11valp1[i] = x[jt + nY*(Time+1) + 2*i*(Time+1) + 1];\n')
f.write('        K11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i)*Time + jt];\n')
f.write('        dK11val[i] = x[jt + (nY+2*i+1)*(Time+1)];\n')
f.write('        dK11valp1[i] = x[jt + (nY+2*i+1)*(Time+1)+1];\n')
f.write('        dK11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i+1)*Time + jt];\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*jt];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdval2[%d] = %s[2*jt+1];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdvalp1[%d] = %s[2*jt+2];\n' % (i,discAzero.Ldata[i]))

for i in range(nI):
    f.write('     Ival[%d] = %s[2*jt];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ival2[%d] = %s[2*jt+1];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ivalp1[%d] = %s[2*jt+2];\n' % (i,discAzero.Lstimuli[i]))

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];\n')
f.write('     } //end for loop\n')
f.write('\n')

VObj = discAzero.VObj
data = len(VObj[0])
model = len(VObj[1])
# Do the gradient elements for the data aspect of the objective function
for i in range(data):
    if VObj[0][i][2] < nY+2*nU:
	f.write('    grad_f[jt+%d*(Time+1)] += (%s)/(2*Time+1);\n' % (VObj[0][i][2], VObj[0][i][0]))
    elif VObj[0][i][2] < 2*(nY+2*nU):
	f.write('    grad_f[jt+1+%d*(Time+1)] += (%s)/(2*Time+1);\n' % (VObj[0][i][2] - (nY+2*nU), VObj[0][i][0]))
    elif VObj[0][i][2] < 3*(nY+2*nU):
	f.write('    grad_f[(Time+1)*(2*nU+nY) + %d*Time + jt] += (%s)/(2*Time+1);\n' % (VObj[0][i][2] - 2*(nY+2*nU), VObj[0][i][0]))
    else:
	f.write('     grad_f[(2*Time+1)*(nY+2*nU)+%d] += (%s)/(2*Time+1);\n' % (VObj[0][i][2]-3*(nY+2*nU), VObj[0][i][0]))
f.write('\n')
#  Now do the model aspect gradient
for i in range(model):
    eqnno = VObj[1][i][2]
    if eqnno < nY+2*nU:
	f.write('    grad_f[jt+%d*(Time+1)] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2], VObj[1][i][1]-1, VObj[1][i][0]))
    elif eqnno < 2*(nY+2*nU):
	f.write('    grad_f[jt+1+%d*(Time+1)] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2]-(nY+2*nU), VObj[1][i][1]-1, VObj[1][i][0]))
    elif eqnno < 3*(nY+2*nU):
	f.write('    grad_f[(Time+1)*(2*nU+nY) + %d*Time + jt] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2]-2*(nY+2*nU), VObj[1][i][1]-1, VObj[1][i][0]))
    else:
	f.write('     grad_f[(2*Time+1)*(nY+2*nU)+%d] += bounds[%s][3]*(%s)/(2*Time+1);\n' % (VObj[1][i][2]-3*(nY+2*nU), VObj[1][i][1]-1,VObj[1][i][0]))
f.write('\n')
f.write('  } //end for loop\n\n')

# Add code for last gradient element

f.write('// Add last element\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[Time + i*(Time+1)];\n')
f.write('        Xvalp1[i] = 0;\n')
f.write('        Xval2[i] = 0;\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index i=0;i<nU;i++) {\n')
f.write('        K11val[i] = x[Time + nY*(Time+1) + 2*i*(Time+1)];\n')
f.write('        K11valp1[i] = 0;\n')
f.write('        K11val2[i] = 0;\n')
f.write('        dK11val[i] = x[Time + (nY+2*i+1)*(Time+1)];\n')
f.write('        dK11valp1[i] = 0;\n')
f.write('        dK11val2[i] = 0;\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*Time];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdval2[%d] = 0;\n' % i)
    f.write('     Xdvalp1[%d] = 0;\n' % i)

for i in range(nI):
    f.write('     Ival[%d] = %s[2*Time];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ival2[%d] = 0;\n' % i)
    f.write('     Ivalp1[%d] = 0;\n' % i)

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')

for i in range(data):
    if VObj[0][i][2] < nY+2*nU:
	f.write('    grad_f[Time+%d*(Time+1)] += (%s)/(2*Time+1);\n' % (VObj[0][i][2], VObj[0][i][0]))

f.write('\n')
f.write('  return true;\n\
}\n\n\n')



# EVAL_G

f.write('// return the value of the constraints: g(x)\n\
bool %s_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)\n\
{\n' % probu)
f.write('  assert(n == %d*Time+%d);\n' % (alpha, beta))
f.write('  assert(m == 0);\n\n')


f.write('  return true;\n\
}\n\n\n')



# EVAL_JAC_G

f.write('// return the structure or values of the jacobian\n\
bool %s_NLP::eval_jac_g(Index n, const Number* x, bool new_x,\n\
                            Index m, Index nele_jac, Index* iRow, Index* jCol,\n\
                            Number* values)\n\
{\n\n\
return true;\n\
}\n\n\n' % probu)



# EVAL_H

f.write('// return the structure or values of the hessian\n\
bool %s_NLP::eval_h(Index n, const Number* x, bool new_x,\n\
                       Number obj_factor, Index m, const Number* lambda,\n\
                       bool new_lambda, Index nele_hess, Index* iRow,\n\
                       Index* jCol, Number* values)\n\
{\n\n\
if (values == NULL) {\n\
   // return the structure.  This is a symmetric matrix, fill in the lower left\n\
   // triangle only.\n\n' % probu)

# Set up dictionaries to mark start point and length of each Hessian entry
# dictlength will be a string of either 1, Time, or Time +1 depending on how
# many entries in the Hessian there are for an element
# dictstart determines where each of the parts start based on what has come before

f.write('   // Each non-zero Hessian element has its own explicit loop\n\
   // since each element needs a different number of matrix elements\n\n')

dictstart = {}
dictlength = {}

VHes = discAzero.VHes

# sma, med, lar track how many 1, T, T+1 there are
sma = 0
med = 0
lar = 0
# count keeps track of the unique index for each entry
# Start at 10 to be sure it is unique
count = 1
for i in range(cat):
  for j in range(i+1):
      # Only take the unique entries: no xp1/xp1 for instance
      if Hesindex[i][j] == -1:
         row = i
	 col = j
	 temp = Hessize[i][j]
	 start = '%d*(Time+1)+%d*(Time)+%d' % (lar,med,sma)
         if temp == 's':
            length = '1'
	    sma = sma + 1
         elif temp == 'm':
            length = 'Time'
	    med = med + 1
         elif temp == 'l':
            length = 'Time+1'
	    lar = lar + 1 
         d1 = {count:start}  
         d2 = {count:length}
	 # Update the Hesindex matrix to reflect count
	 Hesindex[i][j] = count
	 count = count + 1
	 dictstart.update(d1)
	 dictlength.update(d2)
	 f.write('\n   for(Index jt=0;jt<%s;jt++) {\n' % length)
	 f.write('     iRow[%s+jt] = ' % start)
# Now set out which row and column the derivatives are in respect to
         if row < 2*(nY+2*nU):
	    if row % 2 == 0:
		f.write('(Time+1)*%d+jt;\n' % (row/2))
	    else:
		f.write('(Time+1)*%d+jt+1;\n' % (row/2))
	 elif row < 3*(nY+2*nU):
		f.write('(Time+1)*%d+Time*%d+jt;\n' % (nY+2*nU,row-2*(nY+2*nU)))
	 else:
	    f.write('2*Time*%d+%d;\n' % (nY+2*nU, row-2*(nY+2*nU)))
	 f.write('     jCol[%s+jt] = ' % start)
         if col < 2*(nY+2*nU):
	    if col % 2 == 0:
		f.write('(Time+1)*%d+jt;\n' % (col/2))
	    else:
		f.write('(Time+1)*%d+jt+1;\n' % (col/2))
	 elif col < 3*(nY+2*nU):
		f.write('(Time+1)*%d+Time*%d+jt;\n' % (nY+2*nU,col-2*(nY+2*nU)))
	 else:
	    f.write('2*Time*%d+%d;\n' % (nY+2*nU, col-2*(nY+2*nU)))
	 f.write('   }\n')
      else:
         # Change Hesindex to reflect the correct index
	 row = i
	 col = j
	 if i < 2*(nY+2*nU):
	    Hesindex[i][j] = Hesindex[i-1][j-1]
	 else:
	    Hesindex[i][j] = Hesindex[i][j-1]
f.write('}\n\n\
else {\n\
  // return the values.  This is a symmetric matrix, fill the lower left\n\
  // triangle only\n\
  // initialize the values array\n\
  // Point to the initial starting spot for the Hessian elements\n\n\
  for(Index jt=0;jt<%d*(Time+1)+%d*Time+%d;jt++) values[jt] = 0.; // Initialize matrix' % (lar, med ,sma))
f.write('\n\n   // fill the objective portion\n\n')

# Loop over all other elements

f.write('  for(Index jt=0;jt<Time;jt++) {\n\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[jt + i*(Time+1)];\n')
f.write('        Xvalp1[i] = x[jt + i*(Time+1) + 1];\n')
f.write('        Xval2[i] = x[(Time+1)*(nY+2*nU) + jt + i*(Time)];\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index i=0;i<nU;i++) {\n')
f.write('        K11val[i] = x[jt + nY*(Time+1) + 2*i*(Time+1)];\n')
f.write('        K11valp1[i] = x[jt + nY*(Time+1) + 2*i*(Time+1) + 1];\n')
f.write('        K11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i)*Time + jt];\n')
f.write('        dK11val[i] = x[jt + (nY+2*i+1)*(Time+1)];\n')
f.write('        dK11valp1[i] = x[jt + (nY+2*i+1)*(Time+1)+1];\n')
f.write('        dK11val2[i] = x[(Time+1)*(nY+2*nU) + (nY+2*i+1)*Time + jt];\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*jt];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdval2[%d] = %s[2*jt+1];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdvalp1[%d] = %s[2*jt+2];\n' % (i,discAzero.Ldata[i]))

for i in range(nI):
    f.write('     Ival[%d] = %s[2*jt];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ival2[%d] = %s[2*jt+1];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ivalp1[%d] = %s[2*jt+2];\n' % (i,discAzero.Lstimuli[i]))

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];\n')
f.write('     } //end for loop\n')
f.write('\n')
# Do singletons first
for i in range(len(VHes)):
   Hesrow = VHes[i][1]
   Hescol = VHes[i][2]
   if Hessize[Hesrow][Hescol] == 's':
      constraint = VHes[i][0]
      count = Hesindex[Hesrow][Hescol]
      string = VHes[i][4]
      start = dictstart[count]
      if constraint == 0:
          R = '1'
      else:
          R = 'bounds[%s][3]' % (constraint - 1) 
      f.write('    values[%s] += ' % start)
      f.write('obj_factor*%s*(%s)/(2*Time+1);\n\n' % (R,string))
#  Check if this works correctly!
for i in range(len(VHes)):
    Hesrow = VHes[i][1]
    Hescol = VHes[i][2]
    if Hessize[Hesrow][Hescol] != 's':
       constraint = VHes[i][0]
       count = Hesindex[Hesrow][Hescol]
       string = VHes[i][4]
       start = dictstart[count]
       if constraint == 0:
          R = '1'
       else:
          R = 'bounds[%s][3]' % (constraint - 1) 
       toggle = VHes[i][3]
       if toggle == 0:	    
	  f.write('   values[%s+jt] += ' % start)
	  f.write('obj_factor*%s*(%s)/(2*Time+1);\n\n' % (R, string))
       elif toggle == 1:
          newstart = start + '+ 1'
	  f.write('   values[%s+jt] += ' % newstart)
	  f.write('obj_factor*%s*(%s)/(2*Time+1);\n\n' % (R, string))

f.write('   } // end for loop \n\n')

#  Add in code for last element from data part
f.write('// Add last element\n')
f.write('     for(Index i=0;i<nY;i++) {\n')
f.write('        Xval[i] = x[Time + i*(Time+1)];\n')
f.write('        Xvalp1[i] = 0;\n')
f.write('        Xval2[i] = 0;\n')
f.write('     } //end for loop\n\n')

f.write('     for(Index i=0;i<nU;i++) {\n')
f.write('        K11val[i] = x[Time + nY*(Time+1) + 2*i*(Time+1)];\n')
f.write('        K11valp1[i] = 0;\n')
f.write('        K11val2[i] = 0;\n')
f.write('        dK11val[i] = x[Time + (nY+2*i+1)*(Time+1)];\n')
f.write('        dK11valp1[i] = 0;\n')
f.write('        dK11val2[i] = 0;\n')

f.write('     } //end for loop\n\n')

for i in range(nM):
    f.write('     Xdval[%d] = %s[2*Time];\n' % (i,discAzero.Ldata[i]))
    f.write('     Xdval2[%d] = 0;\n' % i)
    f.write('     Xdvalp1[%d] = 0;\n' % i)

for i in range(nI):
    f.write('     Ival[%d] = %s[2*Time];\n' % (i,discAzero.Lstimuli[i]))
    f.write('     Ival2[%d] = 0;\n' % i)
    f.write('     Ivalp1[%d] = 0;\n' % i)

f.write('\n')
f.write('     for(Index i=0;i<nP;i++) {\n')
f.write('        Pval[i] = x[(2*Time+1)*(nY+2*nU)+i];\n')
f.write('     } //end for loop\n\n')
f.write('\n')

for i in range(len(VHes)):
    constraint = VHes[i][0]
    if constraint == 0:
        Hesrow = VHes[i][1]
        Hescol = VHes[i][2]
        if Hessize[Hesrow][Hescol] == 'l':
           string = VHes[i][4]
           count = Hesindex[Hesrow][Hescol]
           start = dictstart[count]
           newstart = start + '+ Time'
           f.write('   values[%s] += ' % newstart)
           f.write('obj_factor*(%s)/(2*Time+1);\n\n' % string)

f.write('  } // end else \n\n')



f.write('   return true;\n\
}\n\n\n')


# FINALIZE_SOLUTION


f.write('\
void %s_NLP::finalize_solution(SolverReturn status,\n\
                        Index n, const Number* x, const Number* z_L, const Number* z_U,\n\
                        Index m, const Number* g, const Number* lambda,\n\
                        Number obj_value,\n\
                        const IpoptData* ip_data,\n\
                        IpoptCalculatedQuantities* ip_cq)\n\
{\n\
  // here is where the solution is written to file\n\n' % probu)


f.write('  FILE *OUTPUT1;\n')
f.write('  FILE *OUTPUT2;\n')
f.write('  FILE *OUTPUT4;\n')
f.write('  FILE *OUTPUT5;\n')

temp1 = "%e"
temp2 = "%d"
temp3 = "\\n"
temp4 = "\t"



f.write('\n\
  OUTPUT1 = fopen ("param.dat","w");\n\
  OUTPUT2 = fopen ("data.dat","w");\n\
  OUTPUT4 = fopen ("Obj.dat","w");\n\n\
  OUTPUT5 = fopen ("carl.input","w");\n\n\
  // Final parameters\n')
for i in range(nP):
   f.write('\
    fprintf (OUTPUT1, "%s%s%s%s%s%s%s%s", bounds[nY+2*nU+%d][0],bounds[nY+2*nU+%d][1],x[(2*Time+1)*(nY+2*nU)+%d]);\n\
    printf("Parameter[%d] = %s%s", x[(2*Time+1)*(nY+2*nU)+%d]);\n\n' % (discAzero.Lparams[i],temp4,temp1,temp4,temp1,temp4,temp1,temp3,i,i,i,i+1,temp1,temp3,i))


f.write(' // Solution of the primal variables, x')
f.write('\n\n\
  for (Index i=0;i<Time;i++) {\n\
    fprintf(OUTPUT2,"%s ' % temp2)
for i in range(nY+nU+nM):
  f.write('%s ' % temp1)
f.write('%s", 2*i, ' % temp3)
for i in range(nY):
  f.write('x[%d*(Time+1)+i], ' % i)
for i in range(nU):
  f.write('x[(nY+2*%d)*(Time+1)+i], ' % i)
for i in range(nM-1):
  f.write('%s[2*i], ' % discAzero.Ldata[i])
f.write('%s[2*i]);\n' % discAzero.Ldata[nM-1])

f.write('\n\
    fprintf(OUTPUT2,"%s ' % temp2)
for i in range(nY+nU+nM):
  f.write('%s ' % temp1)
f.write('%s", 2*i+1, ' % temp3)
for i in range(nY):
  f.write('x[(nY+2*nU)*(Time+1)+%d*Time+i], ' % i)
for i in range(nU):
  f.write('x[(nY+2*nU)*(Time+1)+(nY+2*%d)*Time+i], ' % i)
for i in range(nM-1):
  f.write('%s[2*i+1], ' % discAzero.Ldata[i])
f.write('%s[2*i+1]);\n' % discAzero.Ldata[nM-1])


f.write('  }\n')

# Last time step

f.write('\n\
    fprintf(OUTPUT2,"%s ' % temp2)
for i in range(nY+nU+nM):
  f.write('%s ' % temp1)
f.write('%s", 2*Time, ' % temp3)
for i in range(nY):
  f.write('x[%d*(Time+1)+Time], ' % i)
for i in range(nU):
  f.write('x[(nY+2*%d)*(Time+1)+Time], ' % i)
for i in range(nM-1):
  f.write('%s[2*Time], ' % discAzero.Ldata[i])
f.write('%s[2*Time]);\n\n' % discAzero.Ldata[nM-1])
#####################################
# Add output for use in cudaPIMC
f.write(' // Output for cudaPIMC')
f.write('\n\n\
  for (Index i=0;i<Time;i++) {\n\
     for (Index j=0;j<nY;j++) {\n\
        fprintf(OUTPUT5,"%s%s", x[j*(Time+1)+i]);\n\n' %(temp1,temp3))
f.write('        }\n')
f.write('      for (Index j=0;j<nY;j++) {\n\
        fprintf(OUTPUT5,"%s%s", x[(nY+2*nU)*(Time+1)+j*Time+i]);\n\n' %(temp1,temp3))
f.write('        }\n')
f.write('     }\n')
f.write('  for (Index j=0;j<nY;j++) {\n\
     fprintf(OUTPUT5,"%s%s", x[j*(Time+1)+Time]);\n\n' %(temp1,temp3))
f.write('     }\n')
f.write('  for (Index j=0;j<nP;j++) {\n\
     fprintf(OUTPUT5,"%s%s", x[(2*Time+1)*(nY+2*nU)+j]);\n\n' %(temp1,temp3))
f.write('     }\n')
        
##################################### 
f.write('\n\n')
f.write(' // Output for Objective value')
f.write('  printf("%s%sObjective value%s");\n' % (temp3,temp3,temp3))
f.write('  printf("f(x*) = %s%s", obj_value);\n\n' % (temp1, temp3))
f.write('  fprintf(OUTPUT4,"%s%s", obj_value);\n\n' % (temp1, temp3))


f.write('  fclose (OUTPUT1);\n')
f.write('  fclose (OUTPUT2);\n')
f.write('  fclose (OUTPUT4);\n')
f.write('  fclose (OUTPUT5);\n')

f.write('}\n')       



f.close( )
