#####################################################
#
#  20 October 2009
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
# 
#  This script writes the program file for a C++ IPOPT 
#  program defined by the vector field in the file
#  equations.txt and written by the script makecode.py.
#  This file creates an instance of the non-linear file
#  defined by makehpp.py and makecode.py, and interfaces
#  with the IPOPT libraries necessary for solving the
#  optimization problem.
#
#  This script has been developed as part of a suite of 
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is 
#  generally applicable to any application needing
#  discAzerod derivatives of a vector field.
#
######################################################

import discAzero

prob = discAzero.Problem
probu = prob.upper()
probl = prob.lower()

FILE = probl + 'minAzero_main.cpp'

new = '\\n'
# The name of the IPOPT main file
f = open(FILE,'w')

f.write('// %s_main.cpp\n\
// Main file for use with IPOPT\n' % (probl))

f.write('\
\n\
#include "IpIpoptApplication.hpp"\n\
#include "%sminAzero_nlp.hpp"\n\
\n\
' % probl)
# this is no longer supported
#// for printf\n\
#ifdef HAVE_CSTDIO\n\
# include <cstdio>\n\
#else\n\
# ifdef HAVE_STDIO_H\n\
#  include <stdio.h>\n\
# else\n\
#  error "don\'t have header file for stdio"\n\
# endif\n\
#endif\n\
f.write('\
#include <cstdio>\n\
\n\
using namespace Ipopt;\n\
\n\
int main(int argv, char* argc[])\n\
{\n\
\n\
  // Create a new instance of your nlp\n\
  //  (use a SmartPtr, not raw)\n\
  SmartPtr<TNLP> mynlp = new %s_NLP();\n\
\n\
  // Create a new instance of IpoptApplication\n\
  //  (use a SmartPtr, not raw)\n\
  SmartPtr<IpoptApplication> app = new IpoptApplication();\n\
\n\
  // Change some options\n\
  // Note: The following choices are only examples, they might not be\n\
  //       suitable for your optimization problem.\n\
  app->Options()->SetNumericValue("tol", 1e-12);\n\
  app->Options()->SetStringValue("mu_strategy", "adaptive");\n\
  app->Options()->SetStringValue("output_file", "ipopt.out");\n\
  // The following overwrites the default name (ipopt.opt) of the\n\
  // options file\n\
  app->Options()->SetStringValue("option_file_name", "%s.opt");\n\
\n\
  // Intialize the IpoptApplication and process the options\n\
  ApplicationReturnStatus status;\n\
  status = app->Initialize();\n\
  if (status != Solve_Succeeded) {\n\
    printf("%s%s *** Error during initialization!%s");\n\
    return (int) status;\n\
  }\n\
\n\
  // Ask Ipopt to solve the problem\n\
  status = app->OptimizeTNLP(mynlp);\n\
\n\
  if (status == Solve_Succeeded) {\n\
    printf("%s%s*** The problem solved!%s");\n\
  }\n\
  else {\n\
    printf("%s%s*** The problem FAILED!%s");\n\
  }\n\
\n\
  // As the SmartPtrs go out of scope, the reference count\n\
  // will be decremented and the objects will automatically\n\
  // be deleted.\n\
\n\
  return (int) status;\n\
}\n' % ( probu, probl,new,new,new,new,new,new,new,new,new))
f.close()
