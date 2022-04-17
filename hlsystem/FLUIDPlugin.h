

#ifndef __FLUID_PLUGIN_h__
#define __FLUID_PLUGIN_h__

//#include <GEO/GEO_Point.h>
//
#include <SOP/SOP_Node.h>
#include "fluid_system.h"

namespace HDK_Sample {
class SOP_Fluid : public SOP_Node
{
public:
    static OP_Node		*myConstructor(OP_Network*, const char *,
							    OP_Operator *);

    /// Stores the description of the interface of the SOP in Houdini.
    /// Each parm template refers to a parameter.
    static PRM_Template		 myTemplateList[];

    /// This optional data stores the list of local variables.
    static CH_LocalVariable	 myVariables[];

protected:

    SOP_Fluid(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_Fluid();

    /// Disable parameters according to other parameters.
    virtual unsigned		 disableParms();


    /// cookMySop does the actual work of the SOP computing, in this
    /// case, a LSYSTEM
    virtual OP_ERROR		 cookMySop(OP_Context &context);

    /// This function is used to lookup local variables that you have
    /// defined specific to your SOP.
    virtual bool		 evalVariableValue(
				    fpreal &val,
				    int index,
				    int thread);
    // Add virtual overload that delegates to the super class to avoid
    // shadow warnings.
    virtual bool		 evalVariableValue(
				    UT_String &v,
				    int i,
				    int thread)
				 {
				     return evalVariableValue(v, i, thread);
				 }

private:
    /// The following list of accessors simplify evaluating the parameters
    /// of the SOP.

    // PUT YOUR CODE HERE
	// Here you need to declare functions which need to be called from the .C file to 
	// constantly update the cook function, these functions help you get the current value that the node has
	// Example : To declare a function to fetch angle you need to do it this way 
    exint   CONSTARIANT_ITERATION(exint t) { return evalInt("constraintIteration", 0, t); }
    fpreal  CONSTRAINT_STIFFNESS(fpreal t) { return evalFloat("constraintStiffness", 0, t); }

    fpreal  VISCOSITY(fpreal t) { return evalFloat("viscosity", 0, t); }
    fpreal  VORTICITY_CONFINEMENT(fpreal t) { return evalFloat("vorticityConfinement", 0, t); }
    exint TIME_FRAME(exint t) { return evalInt("timeFrame", 0, t); }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Member variables are stored in the actual SOP, not with the geometry
    /// In this case these are just used to transfer data to the local 
    /// variable callback.
    /// Another use for local data is a cache to store expensive calculations.

	// NOTE : You can declare local variables here like this  
    int		myCurrPoint;
    int		myTotalPoints;
    FluidSystem* myFS;
    std::vector<std::vector<glm::dvec3>> totalPos; // all postions at all frame (0-60)
};
} // End HDK_Sample namespace

#endif
