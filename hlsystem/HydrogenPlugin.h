

#ifndef __HYDROGEN_PLUGIN_h__
#define __HYDROGEN_PLUGIN_h__

//
#include <DOP/DOP_Node.h>
#include "FluidSystem.h"

namespace HDK_Sample {
class DOP_Hydrogen : public DOP_Node
{
public:
    static OP_Node		*myConstructor(OP_Network*, const char *name,
							    OP_Operator *op);

    /// Stores the description of the interface of the SOP in Houdini.
    /// Each parm template refers to a parameter.
    static PRM_Template		 myTemplateList[];

    /// This optional data stores the list of local variables.
    static CH_LocalVariable	 myVariables[];

protected:

    DOP_Hydrogen(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~DOP_Hydrogen();

    ///// Disable parameters according to other parameters.
    //virtual unsigned		 disableParms();


    ///// cookMySop does the actual work of the SOP computing, in this
    ///// case, a LSYSTEM
    //virtual OP_ERROR		 cookMySop(OP_Context &context);

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
    // 
    // // our own parameters
    // particle separation
    // constraint iteration --------
    // constraint stiffness -----------
    // Maximum acceleration
    // Tensile Radius - artificial force kernel radius
    // Viscosity ---------
    // Vorticity Confinement ----------
    // 
    // kernel radius
    // Maximum correction
    // 
    exint   CONSTARIANT_ITERATION(exint t){ return evalInt("constrantIteration", 0, t); }
    fpreal  CONSTRAINT_STIFFNESS(fpreal t) { return evalFloat("constraintStiffness", 0, t); }

    fpreal  VISCOSITY(fpreal t) { return evalFloat("viscosity", 0, t); }
    fpreal  VORTICITY_CONFINEMENT(fpreal t) { return evalFloat("vorticityConfinement", 0, t); }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Member variables are stored in the actual SOP, not with the geometry
    /// In this case these are just used to transfer data to the local 
    /// variable callback.
    /// Another use for local data is a cache to store expensive calculations.

	// NOTE : You can declare local variables here like this  
    int		myCurrPoint;
    int		myTotalPoints;
    FluidSystem myFS;


    //from DOP_GroupAndApply.h
protected:
    void                 processObjectsSubclass(fpreal time,
        int foroutputidx,
        const SIM_ObjectArray& objects,
        DOP_Engine& engine) override;
    void                 getInputInfoSubclass(int inputidx,
        DOP_InOutInfo& info) const override;
    void                 getOutputInfoSubclass(int inputidx,
        DOP_InOutInfo& info) const override;

private:
    void		 GROUP(UT_String& str, fpreal t);
    int			 INPUTINDEX(fpreal t);
};
} // End HDK_Sample namespace

#endif
